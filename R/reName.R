reName <-
function(file_name, proccesing="norm", scale='linear', resolution='25bp', ext='.wig', prefix='C', uid=1) {	
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
