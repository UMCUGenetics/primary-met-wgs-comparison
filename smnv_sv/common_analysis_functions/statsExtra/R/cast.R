#' Cast a long dataframe into a wide dataframe
#'
#' @param df A long form dataframe to cast
#' @param row.var Column name for the values to map to the rows of the output dataframe
#' @param col.var Column name for the values to map to the columns of the output dataframe
#' @param value.var Column name for the values to map each cell in the output dataframe
#' @param fill.value A value to fill cells. Default: `FALSE`
#' @param strings.as.factors Convert strings to factors?
#' @param verbose Show progress messages?
#'
#' @return A wide form dataframe
#' @export
#'
cast <- function(df, row.var, col.var, value.var, fill.value=FALSE, strings.as.factors=T, verbose=F){
   if(F){
      df=fusions
      row.var='sample'
      col.var='fusion_name'
      value.var='dummy_bool'
      fill.value=FALSE
      verbose=T
   }

   df <- df[,c(row.var,col.var,value.var)]

   ## --------------------------------
   if(verbose){ message('Getting unique values...') }
   getUniqValues <- function(v){
      v <- if(is.factor(v)){
         levels(v)
      } else {
         unique(sort(v))
      }
   }

   row_levels <- getUniqValues(df[,row.var])
   col_levels <- getUniqValues(df[,col.var])
   value_levels <- getUniqValues(df[,value.var])

   values_are_strings <- is.factor(value_levels) | is.character(value_levels)

   ## --------------------------------
   if(verbose){ message('Coercing factors to characters...') }
   ## ... to prevent errors
   charToFactor <- function(v){
      if(is.factor(v)){
         v <- as.character(v)
      }
      return(v)
   }

   df[,value.var] <- charToFactor(df[,value.var])
   df[,col.var] <- charToFactor(df[,col.var])

   ## --------------------------------
   ## Make vectors per sample
   values_template <- structure(rep(fill.value, length(col_levels)), names=col_levels)

   if(verbose){ message('Splitting dataframe by `row.var`...') }
   df_split <- split(df[,c(value.var, col.var)], df[,row.var])

   if(verbose){ message('Assigning values per row...') }
   ## Use as.data.frame() for quick `do.call(cbind, df)`
   m <- as.data.frame(lapply(df_split, function(i){
      exist_values <- structure(i[,value.var], names=i[,col.var])
      out <- values_template
      out[names(exist_values)] <- exist_values
      return(out)
   }))
   m <- as.data.frame(t(m))

   ## Convert characters to factors
   if(strings.as.factors & values_are_strings){
      if(verbose){ message('Converting strings to factors...') }
      m <- as.data.frame(lapply(m, function(i){
         factor(i, unique(c(fill.value,value_levels)))
      }))
   }

   colnames(m) <- col_levels
   rownames(m) <- names(df_split)
   m <- m[row_levels,]

   return(m)
}
