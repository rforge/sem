combineModels <- function (...){
	UseMethod("combineModels")
}

combineModels.semmod <- function(..., warn=TRUE){
	model <- rbind(...)
	removeRedundantPaths(model, warn=warn)
	model
}

update.semmod <- function (object, file = "", ...) {
	delete.model.element <- function(delete.text, old.model, 
		type = "path") {
		type <- match.arg(type, c("path", "variable", "coefficient"))
		col.index <- list(path = 1, variable = 1, coefficient = 2)[[type]]
		delete.text <- strip.white(delete.text)
		old.model <- old.model[which(is.na(pmatch(strip.white(old.model[, 
								col.index]), delete.text))), ]
		class(old.model) <- "semmod"
		return(old.model)
	}
	modmat <- scan(file = file, what = list(change = "", var1 = "", 
			var2 = "", var3 = "", var4 = ""), sep = ",", strip.white = TRUE, 
		comment.char = "#", fill = TRUE)
	modmat <- cbind(modmat$change, modmat$var1, modmat$var2, 
		modmat$var3)
	if ("add" %in% modmat[, 1]) {
		addmat <- modmat[which(modmat[, 1] == "add"), 2:4, drop=FALSE]
		if (length(addmat) == 3) 
			addmat <- matrix(addmat, ncol = 3)
		class(addmat) <- "semmod"
		object <- combineModels(object, addmat, warn=FALSE)
	}
	if ("delete" %in% modmat[, 1]) {
		deletemat <- modmat[which(modmat[, 1] == "delete"), 2:3, drop=FALSE]
		if (length(deletemat) == 2) 
			deletemat <- matrix(deletemat, ncol = 2)
		for (i in 1:length(deletemat[, 1])) object <- delete.model.element(deletemat[i, 
					1], object, deletemat[i, 2])
	}
	if ("replace" %in% modmat[, 1]) {
		submat <- modmat[which(modmat[, 1] == "replace"), 2:3, drop=FALSE]
		object[, 1:2] <- mapply(function(x, y) gsub(x, y, object[, 
						1:2]), submat[, 1], submat[, 2])
	}
	removeRedundantPaths(object)
}