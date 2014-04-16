
require(digest)
reName <- function(file_name, proccesing="BEADS", scale='linear', resolution='1bp', ext='.bw', prefix='C', uid=1) {	
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
		#warning('Wrong name format')
	  pos <- regexpr("\\.([[:alnum:]]+)$", file_name)
		return(sprintf("%s_%s_%s_%s%s", ifelse(pos > -1L, substring(file_name, 0, pos-1), file_name), proccesing, scale, resolution, ext))
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
