CheckName <-
function(file_name) {
	if ( grepl("^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)", file_name, perl=T) ) {
		
		ff <- strsplit( substr(file_name, 1, unlist(gregexpr('\\.\\w+$', file_name))-1) , '(\\^|_)', perl=T)[[1]]
		checksum <- substr(digest(sprintf('%s^%s_%s^%s^%s^%s_%s^%s^%s_%s', ff[1], ff[2], ff[3], ff[4], ff[5], ff[6], ff[7], ff[8], ff[9], ff[10])), 1, 2)
		return(substr(ff[11], 2,3) == checksum)	
	
	} else {
		warning('Wrong name format')
		return(NA)	
	}
}
