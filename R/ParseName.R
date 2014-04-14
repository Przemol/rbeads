ParseName <-
function(file_name) {
	if ( grepl("^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)_([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)\\^([A-Za-z0-9\\-\\.]+)", file_name, perl=T) ) {
		
		if( unlist(gregexpr('\\.\\w+$', file_name)) > 0) {
			ff <- strsplit( substr(file_name, 1, unlist(gregexpr('\\.\\w+$', file_name))-1) , '(\\^|_)', perl=T)[[1]]
			ext <- substr(file_name, unlist(gregexpr('\\.\\w+$', file_name)), nchar(file_name))
		} else {
			ff <- strsplit(file_name, '(\\^|_)', perl=T)[[1]]
			ext <- NA
		}
		ff <- strsplit( substr(file_name, 1, unlist(gregexpr('\\.\\w+$', file_name))-1) , '(\\^|_)', perl=T)[[1]]
		df <- data.frame(Factor=I(ff[1]), Antibody=I(ff[2]), ExtractID=I(ff[3]), Crosslinker=I(ff[4]), Strain=I(ff[5]), Stage=I(ff[6]), Processing=I(ff[7]), Resolution=I(ff[8]), Scale=I(ff[9]), ContactExpID=I(ff[10]), UID=I(ff[11]), Extension=I(ext))
		return(df)
	} else {
		warning('Wrong name format')
	}
}
