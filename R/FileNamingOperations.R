# TODO: Add comment
# 
# Author: przemol
###############################################################################
require(digest)

reName <- function(file_name, proccesing="norm", scale='linear', resolution='25bp', ext='.wig', prefix='C', uid=1) {	
	if ( grepl("^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)", file_name, perl=T) ) {
		
		if( unlist(gregexpr('\\.\\w+$', file_name)) > 0) {
			ff <- strsplit( substr(file_name, 1, unlist(gregexpr('\\.\\w+$', file_name))-1) , '(\\^|_)', perl=T)[[1]]
		} else {
			ff <- strsplit(file_name, '(\\^|_)', perl=T)[[1]]
		}
		out <- sprintf('%s^%s_%s^%s^%s^%s_%s^%s^%s_%s', ff[1], ff[2], ff[3], ff[4], ff[5], ff[6], proccesing, scale, resolution, ff[10])
		md5 <- substr(digest(out, algo='md5'), 1, 2)
		num = sprintf("%05d", uid)
		return(sprintf('%s^%s%s%s%s', out, prefix, md5, num, ext))
	} else {
		warning('Wrong name format')
		return(sprintf("%s_%s^%s^%s%s", unlist(strsplit(file_name, "\\."))[1], proccesing, scale, resolution, ext))
	}
}
ParseName <- function(x) {
	rx <- "^([^\\^_]+)\\^([^\\^_]+)_([^\\^_]+)\\^([^\\^_]+)\\^([^\\^_]+)\\^([^\\^_]+)_([^\\^_]+)\\^([^\\^_]+)\\^([^\\^_]+)_([^\\^_]+)\\^([^\\^_\\.]+)(_[^\\.]+)?(\\..+)$"
	if ( grepl(rx, x) ) {
		out <- unlist( regmatches(x, regexec(rx, x)) )[-1] 
		names(out) <- c("Factor", "Antibody", "ExtractID", "Crosslinker", "Strain", "Stage", "Processing", "Scale", "Resolution", "ContactExpID", "UID", "Comments", "Extension")
		return( data.frame(t(out[c(1:11,13,12)]), stringsAsFactors = FALSE) )		
	} else warning('Wrong name format')
}

#ParseName <- function(file_name) {
#	if ( grepl("^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)", file_name, perl=T) ) {
#		
#		if( unlist(gregexpr('\\.\\w+$', file_name)) > 0) {
#			ff <- strsplit( substr(file_name, 1, unlist(gregexpr('\\.\\w+$', file_name))-1) , '(\\^|_)', perl=T)[[1]]
#			ext <- substr(file_name, unlist(gregexpr('\\.\\w+$', file_name)), nchar(file_name))
#		} else {
#			ff <- strsplit(file_name, '(\\^|_)', perl=T)[[1]]
#			ext <- NA
#		}
#		ff <- strsplit( substr(file_name, 1, unlist(gregexpr('\\.\\w+$', file_name))-1) , '(\\^|_)', perl=T)[[1]]
#		df <- data.frame(Factor=I(ff[1]), Antibody=I(ff[2]), ExtractID=I(ff[3]), Crosslinker=I(ff[4]), Strain=I(ff[5]), Stage=I(ff[6]), Processing=I(ff[7]), Resolution=I(ff[8]), Scale=I(ff[9]), ContactExpID=I(ff[10]), UID=I(ff[11]), Extension=I(ext))
#		return(df)
#	} else {
#		warning('Wrong name format')
#	}
#}

CheckName <- function(file_name) {
	if ( grepl("^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)", file_name, perl=T) ) {
		
		ff <- strsplit( substr(file_name, 1, unlist(gregexpr('\\.\\w+$', file_name))-1) , '(\\^|_)', perl=T)[[1]]
		checksum <- substr(digest(sprintf('%s^%s_%s^%s^%s^%s_%s^%s^%s_%s', ff[1], ff[2], ff[3], ff[4], ff[5], ff[6], ff[7], ff[8], ff[9], ff[10])), 1, 2)
		return(substr(ff[11], 2,3) == checksum)	
	
	} else {
		warning('Wrong name format')
		return(NA)	
	}
}
