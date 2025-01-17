#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# R6 class for a Recipe --------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Recipe <- R6Class(
  classname = "recipe",
  public = list(

    # formula the model formula.
    formula = NULL,
    # term_info information about original and predicted variable.
    term_info = NULL,
    # steps the recipe steps.
    steps = list(),
    # template the data.
    template = list(),
    # levels variable factor levels.
    levels = NULL,
    # orig_lvls original factor levels
    orig_lvls = NULL,
    # retained variables that are retained.
    retained = NA,
    # requirements packages required
    requirements = NULL,
    # tr_info training info
    tr_info = NULL,
    # whether the model is trained.
    trained = NULL,

    # bake time for each step
    time_bake = NULL,
    # prep time for each step
    time_prep = NULL,

    template_step = NULL,

    # result list that holds the created model features.
    result = list(),


    # vars The variables available from the provided data set.
    vars = NULL,

    initialize = function(formula = NULL, data = NULL, ...) {

      # specify data used with formula notation
      # if (!is.formula(formula)) {
      #   stop("You must specify a valid formula")
      # }

      if (is.null(data) & is.null(formula)) {
        invisble(self)
      }

      if (!any(class(data) %iin% c("list", "data.frame", "data.table", "tbl"))) {
        stop("data must be a data.frame like object or list")
      }

      self$formula <- formula

      # parse the formula
      vars_list <- get_formula_vars(formula = formula, data = unclass(data))
      roles <-  rep.int(
        x = c("predictor", "outcome"),
        times = vapply(vars_list,
                       FUN = length,
                       FUN.VALUE = numeric(1L)
        )
      )
      self$vars <- unlist(vars_list, use.names = FALSE)

      # add variables as a initial step
      self$add_step(StepAddVars$new(terms = self$vars, role = roles))
      self$steps[[1L]]$set_result(unclass(data)[self$vars])
      self$steps[[1L]]$columns <- self$vars

      self$template_step <- 1L

      self$requirements <- list(
        bake = setNames(object = logical(),
                        nm = character())
      )

      invisible(self)
    },


    # @description
    # Add a step to the recipe.
    # @param step The step to add.
    # @return An updated `Recipe` object with a step addded.
    add_step = function(step) {
      self$steps <- append(self$steps, step)
      invisible(self)
    },


    # @description
    # Do prep operations.
    # @param retain retain the step.
    # @return An updated `Recipe` object.
    prep = function(retain = TRUE, steps = NULL) {

      if (is.null(steps)) {
        steps <- 1:length(self$steps)
      }

      self$time_prep <- c()
      nms <- self$get_step_column_names()

      for (i in steps) {
        start_time <- Sys.time()

        # set the columns that are used for the step
        # currently don't specify as a selector
        self$steps[[i]]$columns <- get_terms_from_info(self$steps[[i]]$terms,
                                                       unlist(nms[1:i]))

        # print(self$steps[[i]]$columns)
        # We need to grab the data from previously created steps
        self$steps[[i]]$prep(self$get_step_columns(self$steps[[i]]$columns))


        end_time <- Sys.time()
        elapsed_time <- end_time - start_time
        self$time_prep <- c(self$time_prep, elapsed_time)

      }

      self$retained <- retain

      invisible(self)
    },

    # @description
    # Create the dataset.
    # @param data The input data to the recipe. If it is not specified it uses
    # the data initially provided to the Recipe.
    # @return An updated `Recipe` object with a result that holds a list of
    # features.
    bake = function(data = NULL, steps = NULL) {


      types <- self$get_step_types()

      types_loop <- seq_along(types)

      if (!is.null(steps)) {
        types_loop <- steps
      } else {
        status <- self$get_step_status() # rerun a step
        types_loop <- types_loop[which(status)]
      }

      if (!is.null(data)) {
        self$steps[[self$template_step]]$set_result(unclass(data)[self$vars])
        self$steps[[self$template_step]]$columns <- self$vars
      }


      for (i in types_loop) {

        start_time <- Sys.time()

        columns <- self$steps[[i]]$columns
        if (is.null(columns)) {
          columns <- names(self$result)[self$template_step]
        }

        # The template from step_add_vars is now in self$steps[[1]]

        switch(

          types[i],
          # "add" = self$steps[[i]]$bake(unclass(self$result)[columns]),
          #
          # "modify" = self$steps[[i]]$bake(unclass(self$result)[columns]),

          # "supervise_add" = self$steps[[i]]$bake(unclass(self$result),
          #                                        self$term_info),
          #
          # # "add_from_template" = self$steps[[i]]$bake(unclass(self$template)[columns]),
          #
          # "supervise_augment" = self$steps[[i]]$bake(unclass(self$result),
          #                                            self$term_info,
          #                                            self$steps),
          # default
          # self$steps[[i]]$bake(unclass(self$template)[columns])
          "model" = self$steps[[i]]$bake(self),
          # default
          self$steps[[i]]$bake(self$steps[[self$template_step]])

        )


        end_time <- Sys.time()
        elapsed_time <- end_time - start_time
        self$time_bake <- c(self$time_bake, elapsed_time)


        self$update_term_info(
          step_name = self$steps[[i]]$step_name,
          step_index = i,
          roles = self$steps[[i]]$role
        )

        self$steps[[i]]$rerun <- FALSE

      }

      invisible(self)
    },

    # @description
    # Reduce the recipe to tabular form. Bake and coerce to the desired output
    # type.
    # @param type The output data type: data.frame, data.table, matrix, tibble,
    # @return tabular output of baked Recipe.
    plate = function(type = "df",
                     steps = NULL,
                     step_dont_rerun = NULL,
                     ...) {

      # prep and bake recipe if it hasn't been done
      if (is.null(steps)) {
        steps <- 1:length(self$steps)
      }

      step_prep_bake <- setdiff(steps, step_dont_rerun)
      return_type(x = self$prep(steps = step_prep_bake)$
                    bake(steps = step_prep_bake)$
                    get_result(steps = steps),
                  type = type, ...)

    },
    # @description
    # get info about steps
    # @param type The output data type: data.frame, data.table, matrix, tibble,
    # @return tabular output of baked Recipe.
    tidy = function(type = "df") {

      info <- list()
      for (i in seq_along(self$steps)) {
        info[[i]] <- self$steps[[i]]$tidy(i)
      }

      collapse::rowbind(info)

    },
    # # short summary of training set.
    # # @description
    # # Get info about the training set.
    # # @return data.frame with limited info on the training set
    # train_info = function(x) {
    #   data.frame(
    #     nrows = length(self$template[[1L]])
    #     # ncomplete = collapse::fsum(!collapse::missing_cases(self$template))
    #   )
    # },
    # @description
    # Update the term info after baking a step.
    # @param source Where did the new terms come from: derived, original
    # @param roles The step role.
    # @param type The type of the step: add, modify, delete, etc.
    # @param step_name The name of the step
    # @param step_index The order the step was added
    # @return updated term_info
    update_term_info = function(source = "derived",
                                roles = "predictor",
                                type = "numeric",
                                step_name,
                                step_index) {


      nms <- names(self$result)
      variable <- setdiff(nms, self$term_info$variable)
      variable_rem <- setdiff(self$term_info$variable, nms)

      n_rem <- length(variable_rem)

      n <- length(variable)

      if (length(roles) == 1L) {
        roles <- rep.int(roles, times = n)
      }

      if (length(source) == 1L) {
        source <- rep.int(source, times = n)
      }

      if (n > 0L) {
        self$term_info$variable <- c(self$term_info$variable, variable)
        self$term_info$roles <- c(
          self$term_info$roles,
          roles
        )
        self$term_info$source <- c(
          self$term_info$source,
          source
        )
        self$term_info$type <- c(
          self$term_info$type,
          rep.int(type, times = n)
        )
        self$term_info$sub_type <- c(
          self$term_info$sub_type,
          rep.int(type, times = n)
        )
        self$term_info$step_index <- c(
          self$term_info$step_index,
          rep.int(step_index, times = n)
        )
        self$term_info$step_name <- c(
          self$term_info$step_name,
          rep.int(step_name, times = n)
        )
      }

      if (n_rem > 0L) {
        wh <- which(self$term_info$variable %iin% variable_rem)
        self$term_info$source[wh] <- "removed"
        self$term_info$step_index[wh] <- step_index
      }

      self
    },
    # @description
    # Get the type of the step.
    # @return character vector for the step types
    check_result_lengths = function() {
      n <- collapse::fnunique(collapse::vlengths(self$results))
      if (n > 1L) {
        warning('Result lengths are not equal')
      }
    },

    # @description
    # Get the indices of previously baked steps.
    # @return integer vector of indices
    is_baked = function() {
      if (is.null(self$term_info)) {
        return(FALSE)
      }
      unique(self$term_info$step_index)
    },
    # @description
    # Get the number of steps.
    # @return integer
    n_steps = function() {
      length(self$steps)
    },
    # @description
    # Remove a step
    # @return recipe with a step(s) removed
    remove_steps = function(step_numbers) {
      self$steps <- self$steps[-step_numbers]
      self
    },
    # @description
    # update a step value
    # @return recipe
    update_steps = function(pars, varying) {

      step_nums <- sort(unique(varying$step_number))

      # loop through each step
      for (i in step_nums) {
        wh <- varying[["step_number"]] == i
        p <- pars[wh]
        n <- varying[['name']][wh]

        if (all(n == "coef")) {
          next
        }

        for (j in which(wh)) {
          self$steps[[i]]$update_step(n[j], p[j])
        }

      }
      self
    },
    # @description
    # update parameters for step_nls
    # @return updated parameters
    update_pars = function(pars, varying) {

      n_new_columns <- self$get_step_number_columns()
      is_coef <- self$get_step_is_coef()

      step_nums <- collapse::funique(varying[["step_number"]], sort = TRUE)

      # loop through each step
      p <- vector("list", length(step_nums))
      for (i in step_nums) {

        if (is_coef[i]) {
          wh <- varying[["step_number"]] == i
          p[[i]] <- pars[wh]
        } else {
          n_new <- n_new_columns[i]
          p[[i]] <- rep.int(1.0, n_new)
        }

      }

      unlist(p, use.names = FALSE)
    },
    get_step_is_coef = function() {
      vapply(self$steps, FUN = function(x) x$is_coef(), FUN.VALUE = logical(1L))
    },
    # @description
    # Get the type of the step.
    # @return character vector for the step types
    get_step_types = function() {
      vapply(self$steps, FUN = function(x) x[["type"]], FUN.VALUE = character(1L))
    },
    # @description
    # Get the type of the step.
    # @return character vector for the step types
    get_step_status = function() {
      vapply(self$steps, FUN = function(x) x[["rerun"]], FUN.VALUE = logical(1L))
    },
    # @description
    # Get the type of the step.
    # @return character vector for the step types
    get_varying = function() {
      n <- length(self$steps)
      nm  <- vector("list", n)
      s   <- nm
      u   <- nm
      l   <- nm
      ind <- nm
      for (i in seq_along(self$steps)) {

        self$steps[[i]]$update_varying()
        v  <- self$steps[[i]][["varying"]]

        if (is.null(v)) {
          next
        }

        nm[[i]] <- v[["name"]]
        s[[i]]  <- v[["start"]]
        u[[i]]  <- v[["upper"]]
        l[[i]]  <- v[["lower"]]
        ind[[i]] <- rep.int(i, length(v[["name"]]))

      }

      list(name  = unlist(nm, use.names = FALSE),
           start = unlist(s,  use.names = FALSE),
           upper = unlist(u,  use.names = FALSE),
           lower = unlist(l,  use.names = FALSE),
           step_number = unlist(ind, use.names = FALSE))

    },
    # @description
    # Get the type of the step.
    # @return character vector for the step types
    get_step_field = function(field_name) {
      lapply(self$steps, "[[", field_name)
    },
    # @description
    # Get the result data
    # @return table of results
    get_step_column_names = function() {

      lapply(self$steps, "[[", "columns")

    },
    # @description
    # Get the number of new columns
    # @return table of results
    get_step_number_columns = function() {

      vapply(self$steps,
             FUN = function(x) x$get_number_columns(),
             FUN.VALUE = numeric(1L))

    },
    # @description
    # start of term and length of terms
    # Get the result data
    # @return table of results
    get_term_index = function(nms) {

      nc   <- lapply(self$steps, "[[", "new_columns")

      inds <- lapply(nc, function(x) {

        if (is.null(x)) {
          return(x)
        }

        wh <- nms %iin% x
        if (length(wh) == 0L) {
          return(NULL)
        }

        c(min(wh) - 1L, max(wh) - min(wh) + 1L)

      })

      inds

    },
    # @description
    # Get the outcome data
    # @return table of results
    get_outcome_variable = function(type = "df", steps = NULL) {

      all_vars <- self$steps[[1L]]$columns
      outcome_vars <- all_vars[which(self$steps[[1L]]$role == "outcome")]

      self$steps[[1L]]$result[outcome_vars]

    },
    # @description
    # Get the result data
    # @return table of results
    get_result = function(type = "df", steps = NULL, ...) {

      if (is.null(steps)) {
        steps <- self$template_step:length(self$steps)
      }

      return_type(unlist(lapply(self$steps[steps], "[[", "result"), recursive = FALSE),
                  type = type, ...)

    },
    # @description
    # Get the type of the step.
    # @return character vector for the step types
    get_step_columns = function(column_names) {

      nms <- names(self$steps[[1L]]$result)

      return(self$steps[[1L]]$result[column_names])

    },
    # @description
    # Get responses separated by step
    # @return response table from regression
    get_response_data = function(type = "df") {

      resp <- self$get_step_data("response_data", type = "df")

      if (is.null(resp)) return(NULL)

      if(is.data.frame(resp)) {
        return(return_type(resp, type = type))
      }

      resp <- collapse::rowbind(resp, use.names = FALSE)
      return_type(resp, type = type)

    },
    # @description
    # Get the predictions separated by step
    # @return predicted values based on regression
    get_predict_data = function(type = "df") {

      # at the moment we don't handle multiple ols runs
      fits <- self$get_step_data("fit")[[1]]

      return_type(fits$decomposition, type = type)

    },
    # @description
    # Get the transfer function for steps
    # @return data from a specific step
    get_transfer_data = function(type = "raw") {

      data <- list()
      for (i in seq_along(self$steps)) {
        tmp <- self$steps[[i]][["fft_result"]]

        if (!is.null(tmp)) {
          data[[i]] <- tmp
        }

      }

      data <- data[!sapply(data, is.null)]

      if (type == "raw") {
        return(data)
      }


      if (type %iin% c("df", "dt")) {
        type_name <- "data.frame"

        if (type == "dt") type_name <- "data.table"

        tf <- collapse::rowbind(
          lapply(data, function(z) {
            collapse::rowbind(lapply(z, function(x) {
              dt <- collapse::qDT(x)
              nms <- names(dt)

              if (ncol(dt) == 3L) {
                dt[, variable := nms[1L]]
                setnames(dt, nms[1L], c("value"))
                return(return_type(dt, type = type))
              }

              d <- collapse::pivot(data = dt,
                                   ids = c("frequency", "id"),
                                   how = "longer")
            }))
          }), use.names = TRUE)


        return(return_type(tf, type = type))

      }
    },
    # @description
    # Get the data from a step by name
    # @return data from a specific step
    get_step_data = function(field_name,
                             type = "raw",
                             additional_columns = NULL) {

      data <- list()

      for (i in seq_along(self$steps)) {

        n_list <- self$steps[[i]][[field_name]]

        if (!is.null(n_list)) {


          if (!is.null(additional_columns)) {
            for (j in seq_along(additional_columns)) {

              to_add <- self$steps[[i]][[additional_columns[j]]]

              if (is.null(to_add)) {
                to_add <- NA
              }

              to_add <- list(to_add)
              names(to_add) <- additional_columns[j]
              n_list <- modifyList(x = n_list, val = to_add)
            }
          }
          n_list <- list(n_list)
          data <- append(data, (n_list))

        }
      }

      data <- data[!sapply(data, is.null)]

      if (length(data) == 0L) {
        warning("There were no steps with the provided field_name.")
        return(NULL)
      }
      if (type == "raw") {
        return(data)
      } else {
        return(
          return_type(
            collapse::rowbind(lapply(data, function(x) {
              collapse::qDT(x)}), use.names = FALSE), type)
        )
      }

    },

    # @description
    # Get the time required for prep and bake
    # @return table of elapsed times
    get_elapsed_times = function() {
      data.frame(step_id = sapply(self$steps, "[[", "id"),
                 time_prep = self$time_prep,
                 time_bake = self$time_bake)
    }



  )
)




