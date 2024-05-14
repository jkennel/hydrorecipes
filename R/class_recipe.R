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

      if (!any(class(data) %in% c("list", "data.frame", "data.table", "tbl"))) {
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
      self$template <- unclass(data)

      # add variables as a step
      self$add_step(StepAddVars$new(self$vars, role = roles))$prep()$bake()
      self$term_info$source[] <- "original"

      self$requirements <- list(
        bake = setNames(object = logical(), nm = character())
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
    prep = function(retain = TRUE) {
      # currently this is run twice for the first step
      self$tr_info <- self$train_info()

      # print(self$term_info)

      for (i in seq_along(self$steps)) {
        self$steps[[i]]$prep(unclass(self$template), self$term_info)
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
    bake = function(data = NULL) {

      types <- self$get_step_types()
      baked <- self$is_baked()

      types_loop <- seq_along(types)


      if (is.null(data)) {
        # remove any previously baked
        if (any(baked)) {
          types_loop <- types_loop[-baked]
        }

      } else {
        self$template <- unclass(data)#[unique(self$vars)]
        self$result <- NULL
      }

      for (i in types_loop) {

        columns <- self$steps[[i]]$columns
        if (is.null(columns)) {
          columns <- names(self$result)[1]
        }

        # modify results
        self$result <- switch(
          types[i],

          "add" = append(
            self$result,
            self$steps[[i]]$bake(unclass(self$result)[columns])),

          "modify" = modifyList(
            self$result,
            self$steps[[i]]$bake(unclass(self$result)[columns])),

          "supervise_add" = append(
            self$result,
            self$steps[[i]]$bake(unclass(self$result), self$term_info)),

          "add_from_template" = append(
            self$result,
            self$steps[[i]]$bake(unclass(self$template)[columns])),

          "supervise_augment" = {self$steps[[i]]$bake(unclass(self$result),
                                                      self$term_info,
                                                      self$steps);
            self$result},

          # default
          {self$steps[[i]]$bake(unclass(self$result)[columns]);
            self$result}
        )

        self$update_term_info(
          step_name = self$steps[[i]]$step_name,
          step_index = i,
          roles = self$steps[[i]]$role
        )
      }

      invisible(self)
    },
    # @description
    # Reduce the recipe to tabular form. Bake and coerce to the desired output
    # type.
    # @param type The output data type: data.frame, data.table, matrix, tibble,
    # @return tabular output of baked Recipe.
    plate = function(type = "df") {
      # prep and bake recipe if it hasn't been done
      # if (length(self$result) == 0) {
      self$prep()$bake()
      # }

      return_type(self$result, type = type)
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
    # short summary of training set.
    # @description
    # Get info about the training set.
    # @return data.frame with limited info on the training set
    train_info = function(x) {
      data.frame(
        nrows = length(self$template[[1L]])
        # ncomplete = collapse::fsum(!collapse::missing_cases(self$template))
      )
    },
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

      if (length(roles) == 1) {
        roles <- rep.int(roles, times = n)
      }

      if (length(source) == 1) {
        source <- rep.int(source, times = n)
      }

      if (n > 0) {
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

      if (n_rem > 0) {
        wh <- which(self$term_info$variable %in% variable_rem)
        self$term_info$source[wh] <- "removed"
        self$term_info$step_index[wh] <- step_index
      }

      self
    },
    # @description
    # Get the type of the step.
    # @return character vector for the step types
    get_step_types = function() {
      vapply(self$steps, FUN = function(x) x$type, FUN.VALUE = character(1))
    },
    # @description
    # Get the type of the step.
    # @return character vector for the step types
    check_result_lengths = function() {
      n <- collapse::fnunique(collapse::vlengths(self$results))
      if (n > 1) {
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
    # Get the indices of previously baked steps.
    # @return integer vector of indices
    get_response_data = function(type = "df") {

      resp <- self$get_step_data("response_data")
      resp <- collapse::rowbind(resp)
      return_type(resp, type = type)

    },
    # @description
    # Get the indices of previously baked steps.
    # @return integer vector of indices
    get_predict_data = function(type = "df") {

      # at the moment we don't handle multiple ols runs
      pred <- self$get_step_data("decomposition")
      pred <- collapse::rowbind(pred)
      return_type(pred, type = type)

    },
    # @description
    # Get the indices of previously baked steps.
    # @return integer vector of indices
    get_step_data = function(field_name) {

      data <- lapply(self$steps, function(x) {
        x[[field_name]]
      })

      data[!sapply(data, is.null)]

    }

  )
)




