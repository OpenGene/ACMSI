## MSI Fragments Analysis V1 - 简易版，拥有基本功能，显示不同目标基因的bp范围
## Author： Wen Jiale

## MSI Fragments Analysis V2 - 2023/12/06
## 增加了panel选择和辅助判断，UI界面美化

## MSI Fragments Analysis V3 - 2023/1/09
## 增加了自动判别功能，增加了重比对功能

### Verbatim copy from seqinr package source 3.0-7
read.abif <- function(filename, max.bytes.in.file = file.info(filename)$size, 
	pied.de.pilote = 1.2, verbose = FALSE){
	#
	# Suppress warnings when reading strings with internal nul character:
	#
	RTC <- function(x, ...) suppressWarnings(rawToChar(x, ...))
	#
	# Define some shortcuts:
	#
	SInt32 <- function(f, ...) readBin(f, what = "integer", signed = TRUE, endian = "big", size = 4, ...)
	SInt16 <- function(f, ...) readBin(f, what = "integer", signed = TRUE, endian = "big", size = 2, ...)
	SInt8 <- function(f, ...) readBin(f, what = "integer", signed = TRUE, endian = "big", size = 1, ...)
	UInt32 <- function(f, ...) readBin(f, what = "integer", signed = FALSE, endian = "big", size = 4, ...)
	UInt16 <- function(f, ...) readBin(f, what = "integer", signed = FALSE, endian = "big", size = 2, ...)
	UInt8 <- function(f, ...) readBin(f, what = "integer", signed = FALSE, endian = "big", size = 1, ...)
	f32 <- function(f, ...) readBin(f, what = "numeric", size = 4, ...)
	f64 <- function(f, ...) readBin(f, what = "numeric", size = 8, ...)
	#
	# Load raw data in memory:
	#
	fc <- file(filename, open = "rb")
	rawdata <- readBin(fc, what = "raw", n = pied.de.pilote*max.bytes.in.file)
	if(verbose) {
		print(paste("number of bytes in file", filename, "is", length(rawdata)))
	}
	close(fc)
	#
	# Make a list to store results:
	#
	res <- list(Header = NULL, Directory = NA, Data = NA)
	#
	# Header section is 128 bytes long, located at a fixed position at the
	# beginning of the file. We essentially need the number of item and dataoffset
	#
	res$Header$abif <- RTC(rawdata[1:4])
	if(res$Header$abif != "ABIF") stop("file not in ABIF format")
	if(verbose) print("OK: File is in ABIF format")

	res$Header$version <- SInt16(rawdata[5:6])
	if(verbose) print(paste("File in ABIF version", res$Header$version/100))
	
	res$Header$DirEntry.name <- rawdata[7:10]
	if(verbose) print(paste("DirEntry name: ", RTC(res$Header$DirEntry.name)))

	res$Header$DirEntry.number <- SInt32(rawdata[11:14])
	if(verbose) print(paste("DirEntry number: ", res$Header$DirEntry.number))
	
	res$Header$DirEntry.elementtype <- SInt16(rawdata[15:16])
	if(verbose) print(paste("DirEntry elementtype: ", res$Header$DirEntry.elementtype))

	res$Header$DirEntry.elementsize <- SInt16(rawdata[17:18])
	if(verbose) print(paste("DirEntry elementsize: ", res$Header$DirEntry.elementsize))

	# This one is important:
	res$Header$numelements <- SInt32(rawdata[19:22])
	if(verbose) print(paste("DirEntry numelements: ", res$Header$numelements))
	
	# This one is important too:
	res$Header$dataoffset <- SInt32(rawdata[27:30])
	if(verbose) print(paste("DirEntry dataoffset: ", res$Header$dataoffset))
	dataoffset <- res$Header$dataoffset + 1 # start position is 1 in R vectors
	
	res$Header$datahandle <- SInt32(rawdata[31:34])
	if(verbose) print(paste("DirEntry datahandle: ", res$Header$datahandle))

	res$Header$unused <- SInt16(rawdata[35:128], n = 47)
	# Should be ingnored and set to zero
	res$Header$unused[1:length(res$Header$unused)] <- 0
	if(verbose) print(paste("DirEntry unused: ", length(res$Header$unused), "2-byte integers"))

	#
	# The directory is located at the offset specified in the header,
	# and consist of an array of directory entries.
	# We scan the directory to put values in a data.frame:
	#
	dirdf <- data.frame(list(name = character(0)))
	dirdf$name <- as.character(dirdf$name) # force to characters
	
	for(i in seq_len(res$Header$numelements)){
		deb <- (i-1)*res$Header$DirEntry.elementsize + dataoffset
		direntry <- rawdata[deb:(deb + res$Header$DirEntry.elementsize)]
		dirdf[i, "name"] <- RTC(direntry[1:4])
		dirdf[i, "tagnumber"] <- SInt32(direntry[5:8])
		dirdf[i, "elementtype"] <- SInt16(direntry[9:10])
		dirdf[i, "elementsize"] <- SInt16(direntry[11:12])
		dirdf[i, "numelements"] <- SInt32(direntry[13:16])
		dirdf[i, "datasize"] <- SInt32(direntry[17:20])
		dirdf[i, "dataoffset"] <- SInt32(direntry[21:24])
	}
	if(verbose){
		print("Element found:")
		print(dirdf$name)
	}
	#
	# Save Directory and make a list to store data:
	#
	res$Directory <- dirdf
	res$Data <- vector("list", nrow(dirdf))
	names(res$Data) <- paste(dirdf$name, dirdf$tagnumber, sep = ".")
	#
	# Data extraction:
	#
	for(i in seq_len(res$Header$numelements)){
		deb <- (i-1)*res$Header$DirEntry.elementsize + dataoffset
		# Short data are stored in dataoffset directly:
		if(dirdf[i, "datasize"] > 4){
		debinraw <- dirdf[i, "dataoffset"] + 1
		} else {
		debinraw <- deb + 20
		}
		elementtype <- dirdf[i, "elementtype"]
		numelements <- dirdf[i, "numelements"]
		elementsize <- dirdf[i, "elementsize"]
		data <- rawdata[debinraw:(debinraw + numelements*elementsize)]
		# unsigned 8 bits integer:
		if(elementtype == 1) res$Data[[i]] <- UInt8(data, n = numelements)
		# char or signed 8 bits integer
		if(elementtype == 2){
			res$Data[[i]] <- tryCatch(RTC(data),finally=paste(rawToChar(data,multiple=TRUE),collapse=""),error=function(er){cat(paste("an error was detected with the following  message:",er," but this error was fixed\n",sep=" "))})
		}
		
		# unsigned 16 bits integer:
		if(elementtype == 3) res$Data[[i]] <- UInt16(data, n = numelements)
		# short:
		if(elementtype == 4) res$Data[[i]] <- SInt16(data, n = numelements)
		# long:
		if(elementtype == 5) res$Data[[i]] <- SInt32(data, n = numelements)
		# float:
		if(elementtype == 7) res$Data[[i]] <- f32(data, n = numelements)
		# double:
		if(elementtype == 8) res$Data[[i]] <- f64(data, n = numelements)
		# date:
		if(elementtype == 10)
			res$Data[[i]] <- list(year = SInt16(data, n = 1), 
			month = UInt8(data[-(1:2)], n = 1), day = UInt8(data[-(1:3)], n = 1))
		# time:
		if(elementtype == 11)
			res$Data[[i]] <- list(hour = UInt8(data, n = 1), 
			minute = UInt8(data[-1], n = 1), second = UInt8(data[-(1:2)], n = 1),
			hsecond = UInt8(data[-(1:3)], n = 1))
		# bool:
		if(elementtype == 13) res$Data[[i]] <- as.logical(UInt8(data))
		# pString:
		if(elementtype == 18){
			n <- SInt8(rawdata[debinraw])
			pString <- RTC(rawdata[(debinraw+1):(debinraw+n)])
			res$Data[[i]] <- pString
		}
		# cString:
		if(elementtype == 19) res$Data[[i]] <- RTC(data[1:(length(data) - 1) ])
		# user:
		if(elementtype >= 1024) res$Data[[i]] <- data
		# legacy:
		if(elementtype %in% c(12)) {
			warning("unimplemented legacy type found in file")
		}
		if(elementtype %in% c(6, 9, 14, 15, 16, 17, 20, 128, 256, 384)) {
			warning("unsupported legacy type found in file")
		}
	}
	return(res)
}
# Imports a .fsa file from Applied Biosystems, using the provided converter
read.fsa <- function(
		file,
		lowess = TRUE,
		lowess.top = 200,
		processed = FALSE,
		meta.extra = NULL,
		quiet = FALSE,
		...
	) {
	
	# Parse ABIF
	fsa <- read.abif(file, ...)
	
	# Scan dimensions
	channelCount <- fsa$Data$`Dye#.1`
	if("SCAN.1" %in% names(fsa$Data)) { scanCount <- fsa$Data$SCAN.1
	} else                            { scanCount <- fsa$Data$Scan.1 ### Legacy
	}
	if(scanCount == 0) warning("Warning: No scan stored in this file\n")
	if(channelCount == 0) warning("Warning: No channel stored in this file\n")
	
	# Extract data
	dyeNames <- character(0)
	dyeWavelengths <- integer(0)
	dyeColors <- character(0)
	
	# DATA indexes
	if(is.na(processed))  { processed <- any(sprintf("DATA.%i", c(9:12, 205:299)) %in% names(fsa$Data))
	}
	if(isTRUE(processed)) { dataTracks <- c(9:12, 205:299)
	} else                { dataTracks <- c(1:4, 105:199)
	}
	
	channelValues <- list()
	for(i in 1:channelCount) {
		# Get raw DATA
		values <- fsa$Data[[ sprintf("DATA.%i", dataTracks[i]) ]]
		
		# Apply lowess to reduce time bias
		if(isTRUE(lowess) && scanCount > 0) {
			x <- values
			x[ as.integer(fsa$Data$OfSc.1) + 1L ] <- lowess.top
			x[ x > lowess.top ] <- lowess.top
			values <- values - lowess(x=1:length(values), y=x)$y
		}
		
		# Update value matrix
		channelValues[[i]] <- values

		# Get dye name
		if(sprintf("DyeN.%i", i) %in% names(fsa$Data))                      dyeNames[i] <- fsa$Data[[ sprintf("DyeN.%i", i) ]]
		if(dyeNames[i] == "" && sprintf("DyeS.%i", i) %in% names(fsa$Data)) dyeNames[i] <- fsa$Data[[ sprintf("DyeS.%i", i) ]]
		dyeNames[i] <- gsub(" ", "", dyeNames[i])

		
		# Get dye wavelength
		if(sprintf("DyeW.%i", i) %in% names(fsa$Data)) {
			dyeWavelengths[i] <- fsa$Data[[ sprintf("DyeW.%i", i) ]]
			dyeColors[i] <- wav2RGB(dyeWavelengths[i])
		} else {
			dyeWavelengths[i] <- NA
			dyeColors[i] <- NA
		}
	}
	
	# Channel size consistency
	channelLengths <- sapply(channelValues, length)
	if(isTRUE(processed)) scanCount <- channelLengths[1]
	if(any(channelLengths != scanCount)) {
		stop("Data length inconsistency: ", paste(channelLengths, collapse=", "))
	}
	# Store channels
	x <- matrix(as.double(NA), ncol=channelCount, nrow=scanCount)
	for(i in 1:channelCount) x[,i] <- channelValues[[i]]
	
	# Automatic colors
	if(all(is.na(dyeColors))) {
		dyeColors <- rainbow(channelCount)
		warning("No dye wavelength found, using arbitrary colors")
	}
	
	# Use dye names
	names(dyeWavelengths) <- dyeNames
	names(dyeColors) <- dyeNames
	colnames(x) <- dyeNames
	
	# Attributes
	attr(x, "lowess") <- isTRUE(lowess)
	attr(x, "wavelengths") <- dyeWavelengths
	attr(x, "colors") <- dyeColors
	
	# Metadata to collect
	collect <- c(
		sample = "SpNm",
		well = "TUBE",
		user = "User",
		machine = "MCHN",
		run.name = "RunN",
		run.date = "RUND",
		run.time = "RUNT",
		runModule.name = "RMdN",
		runModule.version = "RMdV",
		runProtocole.name = "RPrN",
		runProtocole.version = "RPrV",
		meta.extra
	)
	
	# Start collection
	meta <- list()
	meta$file <- normalizePath(file)
	for(metaName in names(collect)) {
		values <- fsa$Data[ grep(sprintf("^%s\\.[0-9]+$", collect[ metaName ]), names(fsa$Data)) ]
		if(length(values) > 0L) {
			if(all(sapply(values, is.atomic))) {
				meta[[ metaName ]] <- unlist(values)
				if(length(values) == 1L) names(meta[[ metaName ]]) <- NULL
			} else {
				meta[[ metaName ]] <- as.matrix(as.data.frame(lapply(values, unlist)))
			}
		}
	}
	
	# Reshape dates
	if("runDate" %in% names(meta) && "runTime" %in% names(meta)) {
		dates <- sprintf("%04i-%02i-%02i %02i:%02i:%02i", meta$runDate["year",], meta$runDate["month",], meta$runDate["day",], meta$runTime["hour",], meta$runTime["minute",], meta$runTime["second",])
		names(dates) <- c("runStart", "runStop", "collectionStart", "collectionStop")[ 1:length(dates) ]
		meta$runDate <- NULL
		meta$runTime <- NULL
		meta <- c(meta, dates)
	}
	
	# Injection time (not portable)
	regex <- "^.+<Token>DC_Injection_Time</Token><Value>(.+?)</Value>.+$"
	if("RMdX.1" %in% names(fsa$Data) && grepl(regex, fsa$Data$RMdX.1)) meta$injectionTime <- sub(regex, "\\1", fsa$Data$RMdX.1)
	
	# Store metadata
	attr(x, "metaData") <- meta
	
	# Off scale values (if any)
	attr(x, "offScale") <- as.integer(fsa$Data$OfSc.1) + 1L
	
	# S3 class
	class(x) <- "fsa"
	
	# Print meta data
	if(!quiet) {
		message("FSA metadata : ", paste(sprintf("%s=\"%s\"", names(meta), sapply(meta, "[", 1L)), collapse=", "))
	}
	return(x)
}
# Converts a light wavelength (in nm) into an RGB color namef
# Adapted from : RL <http://codingmess.blogspot.fr/2009/05/conversion-of-wavelength-in-nanometers.html>
wav2RGB <- function(wav) {
	# Initialize
	SSS <- R <- G <- B <- integer(length(wav))
	
	# Color
	zone <- wav >= 380 & wav < 440
	R[zone] <- - (wav[zone] - 440) / (440 - 350)
	B[zone] <- 1
	
	zone <- wav >= 440 & wav < 490
	G[zone] <- (wav[zone] - 440) / (490 - 440)
	B[zone] <- 1
	
	zone <- wav >= 490 & wav < 510
	G[zone] <- 1
	B[zone] <- - (wav[zone] - 510) / (510 - 490)
	
	zone <- wav >= 510 & wav < 580
	R[zone] <- (wav[zone] - 510) / (580 - 510)
	G[zone] <- 1
	
	zone <- wav >= 580 & wav < 645
	R[zone] <- 1
	G[zone] <- - (wav[zone] - 645) / (645 - 580)
	
	zone <- wav >= 645 & wav <= 780
	R[zone] <- 1
	
	# Intensity correction
	zone <- wav >= 380 & wav < 420
	SSS[zone] <- 0.3 + 0.7*(wav[zone] - 350) / (420 - 350)
	
	zone <- wav >= 420 & wav <= 700
	SSS[zone] <- 1
	
	zone <- wav > 700 & wav <= 780
	SSS[zone] <- 0.3 + 0.7*(780 - wav[zone]) / (780 - 700)
	
	# Color name
	out <- rgb(SSS*R, SSS*G, SSS*B)
	
	return(out)
}

try_align <- function(y, trim, allPeaks, fullLadder, useLadder, rMin) {
	rSquared <- 0
	numLadder <- length(useLadder)
	################# 方法一： 从尾到头找
	## 倒序
	tryCatch({
	if (rSquared < rMin) {
		# diff_sequence <- c(20, 20, 20, 20, 6, 14, 20, 20, 20, 20, 6, 14, 20, 20, 20, 20, 6, 14, 20, 20, 10, 10, 20, 6, 14, 20, 20, 20, 20, 6, 14, 20, 20, 20, 20)
		diff_sequence <- rev(diff(useLadder))
		last_value <- tail(allPeaks, 1)
		## 1bp <-> 1idx
		tryPeaks <- allPeaks[y[allPeaks] >= 100]
		### 防止有噪音和太低的情况
		try1trs <- diff(tail(tryPeaks, 2))/20
		try2trs <- diff(tail(allPeaks, 2))/20
		trs <- max(try1trs, try2trs)

		full_sequence <- numeric(length(diff_sequence) + 1)
		candidate_values <- last_value

		for (i in 1:length(full_sequence)) {
			if (!is.null(last_value) && is.integer(last_value) && length(last_value) > 0) {
				full_sequence[i] <- last_value
			} else {
			######## 空集时候直接跳过找下一项
				diff_sequence[i] = diff_sequence[i] + diff_sequence[i-1]
				last_value <- min(full_sequence[full_sequence > 0])
			}
			## 更新汇率
			if (i > 1) {
				trs <- sum(abs(diff(full_sequence))[abs(diff(full_sequence)) < 500])/sum(head(diff_sequence,i-1))
			}
			## 起始峰特殊情况
			if (i < numLadder-1) {
				last_value_range <- c(last_value - (diff_sequence[i] + diff_sequence[i+1] * 0.5) * trs, last_value - diff_sequence[i] * 0.5 * trs)
				candidate_values <- allPeaks[allPeaks >= last_value_range[1] & allPeaks <= last_value_range[2]]
				last_value <- candidate_values[which.max(y[candidate_values])]
			} else {
				last_value_range <- c(last_value - (1.5 *diff_sequence[i]) * trs, last_value - diff_sequence[i] * 0.5 * trs)
				candidate_values <- allPeaks[allPeaks >= last_value_range[1] & allPeaks <= last_value_range[2]]
				if(length(candidate_values) != 0) {
					last_value <- max(candidate_values[y[candidate_values] > 100])
				} else { last_value <- 0 }
			}
		}
		truePeaks <- rev(full_sequence)
		# Restrict modelization on interest zone
		trueSizes <- fullLadder
		usePeaks <- truePeaks[trueSizes %in% useLadder]
		useSizes <- trueSizes[trueSizes %in% useLadder]
		
		# Linear transformation of indexes in bp
		model <- lm(useSizes ~ usePeaks)

		# Check alignment
		rSquared <- summary(model)$adj.r.squared
	}
	}, error = function(e) {
		rSquared <- 0
	})

	################## 方法二：设上下阈值，从原点开始抛
	### 用阈值选取true peaks
	tryCatch({
	if (rSquared < rMin) {
		truePeaks <- allPeaks[y[allPeaks] <= 5000 & y[allPeaks] >= 100]
		if (length(truePeaks) > length(fullLadder)) {
			# Too many peaks retained (artefacts)
			if (trim == "forward") {
				trueSizes <- fullLadder
				truePeaks <- head(truePeaks, length(fullLadder))
			} else if (trim == "backward") {
				trueSizes <- fullLadder
				truePeaks <- tail(truePeaks, length(fullLadder))
			} else {
				stop("Detected more peaks than described in full size ladder", call. = FALSE)
			}
		} else if (length(truePeaks) < length(fullLadder)) {
			# Not enough peaks retained (premature run ending)
			if (trim == "forward") {
				trueSizes <- head(fullLadder, length(truePeaks))
			} else if (trim == "backward") {
				trueSizes <- tail(fullLadder, length(truePeaks))
			} else {
			}
		} else {
			# Fine
			trueSizes <- fullLadder
		}
		# Restrict modelization on interest zone
		usePeaks <- truePeaks[trueSizes %in% useLadder]
		useSizes <- trueSizes[trueSizes %in% useLadder]

		# Linear transformation of indexes in bp
		model <- lm(useSizes ~ usePeaks)
		
		# Check alignment
		rSquared <- summary(model)$adj.r.squared
	}
	}, error = function(e) {
		rSquared <- 0
	})

	###################### 方法三：找到起始峰，对起始峰右边的进行从大到小排序
	tryCatch({
	if (rSquared < rMin) {
		### 用确定起始峰选取true peaks
		last_above_5000 <- max(which(y[allPeaks] > 5000))
		if (is.na(last_above_5000)) {
			cat("No values greater than 5000 found.")
		} else {
			right_of_last_above_5000 <- y[allPeaks][(last_above_5000 + 1):length(y[allPeaks])]
			sorted_right_values <- sort(right_of_last_above_5000, decreasing = TRUE)
			selected_values <- sorted_right_values[1:numLadder]
		}
		selected_indices <- which(y[allPeaks] %in% selected_values)
		truePeaks <- sort(allPeaks[selected_indices])

		trueSizes <- fullLadder

		# Restrict modelization on interest zone
		usePeaks <- truePeaks[trueSizes %in% useLadder]
		useSizes <- trueSizes[trueSizes %in% useLadder]
		
		# Linear transformation of indexes in bp
		model <- lm(useSizes ~ usePeaks)

		# Check alignment
		rSquared <- summary(model)$adj.r.squared
		if (is.na(rSquared)) {
			stop("Unable to compute R-squared (not enough points for a linear model ?)", call. = FALSE)
		}
	}
	}, error = function(e) {
		rSquared <- 0
	})

	##################### 方法四： 找到起始峰，对起始峰右边去除距离相互最近的
	tryCatch({
	if (rSquared < rMin) {
		### 用确定起始峰选取true peaks
		last_above_5000 <- max(which(y[allPeaks] > 5000))
		if (is.na(last_above_5000)) {
			warning("No values greater than 5000 found.")
		} else {
			right_of_last_above_5000 <- y[allPeaks][(last_above_5000 + 1):length(y[allPeaks])]
			selected_right_values <- right_of_last_above_5000[right_of_last_above_5000 > 100]
			selected_values <- sort(selected_right_values, decreasing = TRUE)
			selected_indices <- which(y[allPeaks] %in% selected_values)
			truePeaks <- sort(allPeaks[selected_indices])
			while (length(truePeaks) > length(useLadder)) {
				differences <- abs(diff(truePeaks))
				min_diff_index <- which.min(differences)

				# 找到最小差值是由哪两个相邻的值相减而成的
				element1 <- truePeaks[min_diff_index]
				element2 <- truePeaks[min_diff_index + 1]

				# 比较这两个对应值的y[peaks]
				if (y[element1] > y[element2]) {
					truePeaks <- truePeaks[-(min_diff_index + 1)]
				} else {
					truePeaks <- truePeaks[-min_diff_index]
				}
			}
			# Restrict modelization on interest zone
			usePeaks <- truePeaks[trueSizes %in% useLadder]
			useSizes <- trueSizes[trueSizes %in% useLadder]
		}

	}
	}, error = function(e) {
		rSquared <- 0
	})
	return(list(usePeaks, useSizes))
}
# Size estimation based on few markers searched in predefined ranges
align.fsa <- function(
	x,
	channel = "LIZ",
	fullLadder = c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600),
	useLadder = c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600),
	outThreshold = 0.15,
	noiseLevel = 10,
	surePeaks = 5,
	leakingRatios = c(-1, 10),
	trim = c("forward", "backward", "none"),
	maskOffScale = FALSE,
	rMin = 0.999,
	rescue = FALSE,
	ylim = NULL,
	...
	) {
	# Args
	trim <- match.arg(trim)
	if (is.null(useLadder)) useLadder <- fullLadder

	# Unrecoverable error
	if (!channel %in% colnames(x)) {
		warning("channel \"", channel, "\" not found (available: ", paste(colnames(x), collapse = ", "), ")")
		return("stdNotFound")
	}
	# Postpone errors
	rSquared <- 0
	status <- try(silent = TRUE, expr = {
		# Smoothed channel (offScale masked)
		y <- x[, channel]
		if (isTRUE(maskOffScale)) y[attr(x, "offScale")] <- NA
		y <- ksmooth(x = 1:length(y), y = y, kernel = "normal", bandwidth = 5, n.points = length(y))$y

		# Local maxima (above noise level)
		allPeaks <- which(head(tail(y, -1), -1) >= tail(y, -2) & head(y, -2) < head(tail(y, -1), -1) & head(tail(y, -1), -1) > noiseLevel) + 1L

		# Exclude outliers comparing intensity with sure peaks ('surePeaks' last ones)
		surePeaks.i <- tail(allPeaks, surePeaks)
		sureIntensity <- median(y[surePeaks.i])
		if (outThreshold < 1) {
			threshold <- outThreshold * sureIntensity
		} else {
			threshold <- outThreshold
		}
		outlierPeaks <- allPeaks[abs(y[allPeaks] - sureIntensity) >= threshold]
		truePeaks <- allPeaks[abs(y[allPeaks] - sureIntensity) < threshold]
		## 第一次尝试不清洗数据
		tryCatch({
			align_result <- try_align(y, trim, allPeaks, fullLadder, useLadder, rMin)
			usePeaks <- align_result[[1]] 
			useSizes <- align_result[[2]] 
			# Linear transformation of indexes in bp
			model <- lm(useSizes ~ usePeaks)
			rSquared <- summary(model)$adj.r.squared
		}, error = function(e) {
			rSquared <- 0
		})
		## 第二次尝试清洗数据
		tryCatch({
		if (rSquared < rMin) {
			# Exclude peaks with obvious leakage of one channel into others
			leaking.comp <- apply(x[allPeaks, -grep(channel, colnames(x))] < leakingRatios[1] * x[allPeaks, channel], 1, any)
			leaking.high <- apply(x[allPeaks, -grep(channel, colnames(x))] > leakingRatios[2] * x[allPeaks, channel], 1, any)
			## 尝试清洗comb
			tryPeaks <- allPeaks[!leaking.comp]
			align_result <- try_align(y, trim, tryPeaks, fullLadder, useLadder, rMin)
			usePeaks <- align_result[[1]] 
			useSizes <- align_result[[2]] 
			# Linear transformation of indexes in bp
			model <- lm(useSizes ~ usePeaks)
			rSquared <- summary(model)$adj.r.squared
		}
		}, error = function(e) {
			rSquared <- 0
		})
		tryCatch({
		if (rSquared < rMin) {
				## 尝试清洗high
				tryPeaks <- allPeaks[!leaking.high]
				align_result <- try_align(y, trim, tryPeaks, fullLadder, useLadder, rMin)
				usePeaks <- align_result[[1]] 
				useSizes <- align_result[[2]] 
				# Linear transformation of indexes in bp
				model <- lm(useSizes ~ usePeaks)
				rSquared <- summary(model)$adj.r.squared
		}
		}, error = function(e) {
			rSquared <- 0
		})
		tryCatch({
		if (rSquared < rMin) {
				## 尝试清洗comb和high
				tryPeaks <- allPeaks[!leaking.comp & !leaking.high]
				align_result <- try_align(y, trim, tryPeaks, fullLadder, useLadder, rMin)
				usePeaks <- align_result[[1]] 
				useSizes <- align_result[[2]] 
				# Linear transformation of indexes in bp
				model <- lm(useSizes ~ usePeaks)
				rSquared <- summary(model)$adj.r.squared
		}
		}, error = function(e) {
			rSquared <- 0
		})
	# Store model
	attr(x, "ladderModel") <- coefficients(model)
	attr(x, "ladderExact") <- usePeaks
	attr(x, "ladderR2") <- rSquared
	names(attr(x, "ladderExact")) <- useSizes
	})

	# Alignment rescue data
	if(isTRUE(rescue)) {
		# Replace raw intensities with smoothed ones (for rescue plot only)
		object <- x
		object[,channel] <- y
		
		# Plot profile
		if(is.null(ylim)) ylim <- c(min(y, sureIntensity-threshold, na.rm=TRUE), max(y, sureIntensity+threshold, na.rm=TRUE))
		plot(object, units=ifelse(is(status, "try-error"), "index", "bp"), ladder=!is(status, "try-error"), channels=channel, chanColors="#000000", ylim=ylim, nticks=10, all.bp=FALSE, ...)
		
		# Highlight 'sure peaks'
		if(is(status, "try-error")) { xcor <- surePeaks.i
		} else                      { xcor <- attr(object, "ladderModel")[2] * surePeaks.i + attr(object, "ladderModel")[1]
		}
		points(x=xcor, y=y[surePeaks.i])
		
		# All detected peaks and their status
		if(is(status, "try-error")) { xcor <- allPeaks
		} else                      { xcor <- attr(object, "ladderModel")[2] * allPeaks + attr(object, "ladderModel")[1]
		}
		col <- rep("#888888", length(allPeaks))
		col[ allPeaks %in% outlierPeaks ] <- "#BB3333"
		col[ allPeaks %in% leakingPeaks ] <- "#FFCC00"
		col[ allPeaks %in% truePeaks ] <- "#33BB33"
		if(length(xcor) > 0L) mtext(side=3, at=xcor, text="|", col=col)
		
		# Intensity filter
		xlim <- par()$usr[1:2]
		segments(x0=xlim[1], x1=xlim[2], y0=sureIntensity, y1=sureIntensity, col="#33BB33")
		rect(xleft=xlim[1], xright=xlim[2], ybottom=sureIntensity-threshold, ytop=sureIntensity+threshold, col="#33BB3333", border=NA)
		rect(xleft=xlim[1], xright=xlim[2], ybottom=par("usr")[3], ytop=noiseLevel, col="#BB333333", border=NA)
		
		# Legend
		legend(x="topright", bg="#FFFFFF",
			pch =    c(1,         NA,        NA,          NA,          124,       124,       124,       124,       NA, NA,        NA),
			col =    c("#000000", "#33BB33", NA,          NA,          "#BB3333", "#FFCC00", "#888888", "#33BB33", NA, "#000000", NA),
			lty =    c(NA,        "solid",   NA,          NA,          NA,        NA,        NA,        NA,        NA, "dotted",  NA),
			fill =   c(NA,        NA,        "#33BB3333", "#BB333333", NA,        NA,        NA,        NA,        NA, NA,        NA),
			border = c(NA,        NA,        "#000000",   "#000000",   NA,        NA,        NA,        NA,        NA, NA,        NA),
			legend = c(
				sprintf("'Sure' ladder peaks (%i last)", surePeaks),
				"Ladder intensity (from 'sure' peaks)",
				ifelse(
					outThreshold < 1L,
					sprintf("Tolerance (+/- %g%%)", signif(outThreshold*100L, 3)),
					sprintf("Tolerance (+/- %g)", signif(outThreshold, 3))
				),
				sprintf("Noise (< %g)", signif(noiseLevel, 3)),
				"Excluded : out of tolerance",
				"Excluded : channel leakage",
				"Excluded : trimmed",
				"Retained",
				sprintf("Matching to ladder sizes : %s", trim),
				"Retained for alignment model",
				sprintf("R-squared = %g (%s)", round(rSquared, 6), ifelse(rSquared > rMin, "OK", "MISALIGNED"))
			)
		)
	}
	
	# Call postponed errors
	# if(is(status, "try-error")) {
	# 	stop(conditionMessage(attr(status, "condition")), call.=FALSE)
	# }
	if (any(is.na(usePeaks))) {
		mapping_info = c(attr(x, "metaData")$sample,0)
	} else {
		mapping_info = c(attr(x, "metaData")$sample,round(rSquared, 6))
	}
	return(list(x, mapping_info))
}
# Plot method for "fsa" S3 class
plot.fsa <- function(
		x,
		units = NA,
		channels = NA,
		chanColors = NA,
		ladder = TRUE,
		offScaleCol = "#FF0000",
		offScalePch = "+",
		offScaleCex = 0.4,
		bg = "white",
		fg = "black",
		title = "",
		title.adj = 0,
		title.line = NA,
		xlab = NA,
		ylab = "Intensity",
		xlim = NA,
		ylim = NA,
		xaxt = "s",
		yaxt = "s",
		bty = "o",
		xaxp = NA,
		nticks = 5L,
		all.bp = TRUE,
		peaks.alpha = 48L,
		peaks.srt = 30,
		peaks.adj = c(0, 0),
		peaks.cex = 1.3,
		peaks.font = 2,
		legend.x = "topleft",
		outfile = "channel_data.txt",
		...
	) {
	# Defaults
	if(identical(channels, NA))   channels <- colnames(x)
	if(identical(chanColors, NA)) chanColors <- attr(x, "colors")
	if(identical(units, NA)) {
		if(is.null(attr(x, "ladderModel"))) { units <- "index"
		} else                              { units <- "bp"
		}
	}
	if(identical(xlab, NA))       xlab <- units
	
	# Checks
	if(length(channels) == 0) stop("No channel selected for plotting")
	
	# 'channels' must be a character vector
	if(is.numeric(channels)) channels <- colnames(x)[ channels ]
	
	# Missing channel(s)
	missing <- setdiff(channels, colnames(x))
	if(length(missing) > 0) stop("channel(s) ", paste(missing, collapse=", "), " not found (available: ", paste(colnames(x), collapse=", "), ")")
	
	# 'chanColors' must be a named character vector
	if(is.null(names(chanColors))) {
		# Unamed, supposed to be parallel with channels
		if(length(channels) != length(chanColors)) stop("'chanColors' must be named, or have the same length as 'channels'")
		names(chanColors) <- channels
	} else {
		# Already named, check if all are described
		if(any(! channels %in% names(chanColors))) stop("Some 'channels' elements can't be found in 'chanColors' names")
	}
	
	# Empty plot
	if(nrow(x) == 0 && (identical(xlim, NA) || identical(ylim, NA))) stop("Cannot plot an empty object without 'xlim' and 'ylim'")
	
	# X scale
	units <- match.arg(units, choices=c("index", "bp"))
	if(units == "index") {
		if(nrow(x) > 0) { xcor <- 1:nrow(x)
		} else          { xcor <- integer(0)
		}
		if(identical(xlim, NA)) xlim <- c(1, nrow(x))
		if(identical(xaxp, NA)) xaxp <- NULL
	} else if(units == "bp") {
		if(is.null(attr(x, "ladderModel"))) stop("Can't plot an unaligned object in 'bp' units")
		if(identical(xlim, NA)) xlim <- attr(x, "ladderModel")[2] * c(1, nrow(x)) + attr(x, "ladderModel")[1]
		if(identical(xaxp, NA)) xaxp <- c((xlim[1] %/% nticks)*nticks, (xlim[2] %/% nticks)*nticks, diff(xlim %/% nticks))
		if(nrow(x) > 0) { xcor <- attr(x, "ladderModel")[2] * 1:nrow(x) + attr(x, "ladderModel")[1]
		} else          { xcor <- integer(0)
		}
	}
	
	# ylim
	if(identical(ylim, NA)) {
		# Values to consider
		tmp <- as.vector(x[ xcor >= xlim[1] & xcor <= xlim[2] , channels ])
		tmp <- tmp[ !is.na(tmp) ]
		
		# Max in any
		if(length(tmp) > 0) { ylim <- c(0, max(tmp))
		} else              { ylim <- c(-1, 1)
		}
	}
	
	# Graphical parameters
	savePar <- par(bg=bg, fg=fg, col=fg, col.axis=fg, col.lab=fg, col.main=fg, col.sub=fg)
	on.exit(par(savePar))
	
	# Mask off-scale values
	offScaleValues <- x[ attr(x, "offScale") , , drop=FALSE ]
	x[ attr(x, "offScale") ,] <- NA
	
	# Background
	plot(
		x = NA, y = NA,
		xlim = xlim, ylim = ylim,
		xlab = xlab,
		ylab = ylab,
		xaxt = xaxt,
		yaxt = yaxt,
		bty = bty,
		xaxp = xaxp,
		...
	)
	
	# Peaks
	peaks <- attr(x, "peaks")
	# print(peaks)
	if(!is.null(peaks)) {
		# Ignore invisible peaks
		peaks <- peaks[ !is.na(peaks$color) ,]
		
		# Transparent version of peak colors
		peaks$background <- sprintf(
			"#%s",
			apply(
				rbind(
					col2rgb(peaks$color),
					peaks.alpha
				),
				2,
				function(x) {
					paste(sprintf("%02x", x), collapse="")
				}
			)
		)
		
		# Full height rectangles
		rect(
			xleft = peaks$size.min,
			xright = peaks$size.max,
			ybottom = -1e6,
			ytop = 1e6,
			col = peaks$background,
			border = NA
		)
		
		if(all(c("N0", "N1", "N2") %in% colnames(peaks))) {
			# Allele rectangles
			rect(
				xleft = peaks$size.min,
				xright = peaks$size.max,
				ybottom = c(peaks$N0, peaks$N0) * peaks$height / peaks$normalized,
				ytop = c(peaks$N1, peaks$N2) * peaks$height / peaks$normalized,
				col = peaks$background,
				border = NA
			)
		}
		
		# Names
		text(
			x = (peaks$size.min + peaks$size.max) / 2,
			y = par("usr")[4] + diff(par("usr")[3:4]) / 50,
			labels = rownames(peaks),
			srt = peaks.srt,
			adj = peaks.adj,
			cex = peaks.cex,
			font = peaks.font,
			xpd = NA,
			col = peaks$color
		)
		
		# Genotyping N1-based
		if("present" %in% colnames(peaks)) {
			text(
				x = (peaks$size.min + peaks$size.max) / 2,
				y = par("usr")[4],
				adj = c(0.5, 1),
				col = chanColors[ peaks$channel ],
				labels = ifelse(peaks$present, "+", "-"),
				font = 2
			)
		}
	}

	# Plot channels
	for(h in channels) {
		points(
			x = xcor,
			y = x[,h],
			col = chanColors[h],
			type = "l"
		)

		# 输出每个channel的bp坐标信息和y坐标轴
        cat("Channel:", h, "\n")
        cat("BP Coordinate\tY Coordinate\n")
        for (i in 1:length(xcor)) {
            if (!is.na(x[i, h]) && x[i, h] > 0 && xcor[i] > 0) {
                cat(fsafile, "\t" ,h, "\t", xcor[i], "\t", x[i, h], "\n", file = outfile, append = TRUE)
            }
        }
    }

	scatter_plots <- list()
	for (h in channels) {
		scatter_plot <- plot_ly(
			x = xcor,
			y = x[, h],
			type = "scatter",
			mode = "lines",
			line = list(color = chanColors[h]),
			name = h
			)
			scatter_plots <- c(scatter_plots, list(scatter_plot))
	}

	# Create layout
	layout <- list(
		xaxis = list(range = input$x_range),
		yaxis = list(range = input$y_range),
		dragmode = "zoom",
		title = title,
		showlegend = TRUE
		# Add other layout parameters as needed
	)

	# Combine scatter plots and layout into a plotly object
	plotly_obj <- plot_ly(scatter_plots) %>% layout(layout)

	# Print or save the plotly object as needed
	print(plotly_obj)


	# # Plot bp axis
	# if(units == "bp" && isTRUE(all.bp)) axis(side=1, at=xlim[1]:xlim[2], labels=FALSE)
	
	# Off scale points
	if(length(attr(x, "offScale")) > 0) {
		# Coordinates
		xoff <- xcor[ attr(x, "offScale") ]
		yoff <- offScaleValues
		
		if(length(xoff) > 0) {
			for(h in channels) {
				# Color
				if(is.na(offScaleCol)) { koff <- chanColors[h]
				} else                 { koff <- offScaleCol
				}
				
				# Points
				points(x=xoff, y=yoff[,h], col=koff, type="p", pch=offScalePch, cex=offScaleCex)
			}
		}
	}
	
	# Exact ladder
	if(isTRUE(ladder)) {
		if(!is.null(attr(x, "ladderExact"))) {
			# Add ladder sizes
			at <- switch(units,
				"index" = attr(x, "ladderExact"),
				"bp" = attr(x, "ladderModel")[2] * attr(x, "ladderExact") + attr(x, "ladderModel")[1]
			)
			segments(x0=at, y0=par("usr")[4], x1=at, y1=par("usr")[3]-diff(par("usr")[3:4])/8, lty="dotted", xpd=NA)
			mtext(side=1, at=at, text=names(attr(x, "ladderExact")), line=2)
		} else warning("Can't add ladder without alignment ('ladderExact' attribute)")
	}
	
	# Legend
	inset <- c(par("din")[2]/1000, par("din")[1]/1000)
	legend(legend.x, inset=inset, legend=colnames(x[, channels, drop=FALSE]), col=chanColors[channels], lty="solid", bg="#FFFFFF")
	
	# Title
	title(main=title, adj=title.adj, line=title.line)

	invisible(TRUE)
}
# 添加范围标签的函数
add_range_annotation <- function(p, x_range, y, label, Ccolor) {
	p <- add_trace(
		p,
		x = c(x_range[1], x_range[2]),
		y = y,
		type = "scatter",
		mode = "lines",
		line = list(color = Ccolor, width = 2),
		showlegend = FALSE
	)
	p <- add_annotations(
		p,
		x =  mean(x_range),
		y = y-1000,
		text = label,
		showarrow = FALSE,
		font = list(size = 10)
	)
	return(p)
}

fsa_convert <- function(std_input, files_input, anno_df, session) {
	## 读取输入标准品信息
	std_df <- read.table("www/std_info.txt", header = TRUE, sep = "\t")
	std_select = std_df[std_df$name == std_input, ]
	ladder_str <- strsplit(std_select$ladder, ",")[[1]]
	Inladder <- as.numeric(ladder_str)
	Inchannel <- std_select$channel
	shinyjs::disable("submit")
	peak_data <- data.frame()
	mapping_list = list()

	if (file.exists("www/peak.txt")) {
		file.remove("www/peak.txt")
	}
	fsa_list <- lapply(files_input$datapath, function(file_path) {
		fsa <- read.fsa(file_path)
		channels = colnames(fsa)
		align_result <- align.fsa(fsa, channel = Inchannel, fullLadder = Inladder, useLadder = Inladder, trim = "backward")
		if (length(align_result) != 2) {
			return("stdNotFound")
		}
		fsa <- align_result[[1]]
		mapping_info <- align_result[2]
		mapping_list <<- append(mapping_list, mapping_info)

		xcor <- attr(fsa, "ladderModel")[2] * 1:nrow(fsa) + attr(fsa, "ladderModel")[1]
		chanColors <- attr(fsa, "colors")

		temp_data <- data.frame()

		for(h in channels) {
		temp <- data.frame(sample = attr(fsa, "metaData")$sample,
							channel = h,
							xcor = xcor,
							fsa_value = fsa[, h])
		## 筛选peaks
		temp <- temp[temp$xcor > 30 & temp$fsa_value > 0, ]
		temp_data <- rbind(temp_data, temp)
		}
		peak_data <<- rbind(peak_data, temp_data)
		return(fsa)
	})
	if ("stdNotFound" %in% fsa_list) {
		return("stdNotFound")
	}
	fwrite(peak_data, file = "www/peak.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
	## 比对质量
	mapping_df <- data.frame(matrix(unlist(mapping_list), ncol = 2, byrow = TRUE))
	colnames(mapping_df) <- c("Sample", "MQ")
	write.csv(mapping_df, file = "www/mapping_info.csv", row.names = FALSE)
	shinyjs::enable("submit")
	return(fsa_list)
} 

fsa_plot <- function(fsa_data, selected_colors, new_info, anno_df, input, output, session) {
	## 统一所有样本通道颜色
	color_mapping <- c("blue" = "#0066CC", "green" = "#00CC00", "yellow" = "#FF9900")
	chanColors <- attr(fsa_data[[1]], "colors")
	color_list <- anno_df$Color[match(colnames(fsa_data[[1]]), anno_df$Channel)]
	color_to_channel <- setNames(color_list, colnames(fsa_data[[1]]))
	input_channel <- color_to_channel[match(selected_colors, color_list)]
	input_color <- color_mapping[input_channel]
	input_color <- setNames(input_color, names(input_channel))
	lapply(seq_along(fsa_data), function(i) {
		sampleName = new_info$Sample[i]
		if ((sampleName == '') || is.na(sampleName)) {
			sampleTitle <- attr(fsa_data[[i]], "metaData")$sample
		} else {
			sampleTitle <- paste(sampleName, "-", attr(fsa_data[[i]], "metaData")$sample)
		}
		output[[paste0("plot", i)]] <- renderPlotly({
			fsa <- fsa_data[[i]]
			xcor <- attr(fsa, "ladderModel")[2] * 1:nrow(fsa) + attr(fsa, "ladderModel")[1]
			p <- plot_ly()
			for (h in names(input_channel)) {
				p <- add_trace(
					p,
					x = xcor,
					y = fsa[, h],
					line = list(color = input_color[h]),
					type = "scatter",
					mode = "lines",
					name = h
				)
			}

			# 调用函数添加范围标签
			if (!is.null(input$annofile)) {
				isolate({
					info <- read.table(input$annofile$datapath, header=TRUE, sep=",")
					unique_colors <- unique(info$color)
					height_list <- c(50000)
					for (i in 2:length(unique_colors)) {
						height <- 50000 - (i - 1) * 3000
						height_list <- c(height_list, height)
					}
					for (row in 1:nrow(info)) {
						site <- info$site[row]
						marker <- info$marker[row]
						color <- info$color[row]
						range_str <- info$range[row]
						range <- as.numeric(unlist(strsplit(range_str, "-")))
						height <- height_list[match(color, unique_colors)]
						p <- add_range_annotation(p, range, height, paste(range_str, "bp-", site), color)
					}
				})
			}
			p <- layout(
				p,
				title = sampleTitle,
				xaxis = list(title = "Size", range = input$x_range),
				yaxis = list(title = "RFU", range = input$y_range),
				dragmode = "pan"
			)
			return(p)
		})
	})
}
## 重比对绘制
ladder_plot <- function(fsa_data, new_info, anno_df, input, output, session) {
    lapply(seq_len(length(fsa_data) * 2), function(i) {
		sampleName = new_info$Sample[i/2]
		if (i %% 2 == 0) {
			fsa <- fsa_data[[i/2]]
			## 绘图
			output[[paste0("realign", i)]] <- renderPlotly({
				xcor <- 1:nrow(fsa)
				shapes <- list()
				annotations <- list()
				if ((sampleName == '') || is.na(sampleName)) {
					sampleTitle <- attr(fsa, "metaData")$sample
				} else {
					sampleTitle <- paste(sampleName, "-", attr(fsa, "metaData")$sample)
				}
				for (j in seq_along(attr(fsa, "ladderExact"))) {
					ladder <- attr(fsa, "ladderExact")[[j]]
					name <- names(attr(fsa, "ladderExact"))[j]
					# ladder线
					shape <- list(
					type = "line",
					x0 = ladder, x1 = ladder,
					y0 = 0, y1 = 10000,
					line = list(color = "red", width = 1, dash="dash"),
					text = name
					)
					shapes[[j]] <- shape
					# ladder标签
					annotation <- list( 
					x = ladder,
					y = 10000,
					text = name,
					showarrow = FALSE,
					yshift = 10
					)
					annotations[[j]] <- annotation
				}
				p <- plot_ly(source = paste0("align", i)) %>%
					add_trace(
					x = xcor,
					y = fsa[, 'LIZ'],
					type = "scatter",
					mode = "lines"
					)  %>%
					layout(
						title = sampleTitle,
						shapes = shapes,
						annotations = annotations,
						uirevision = TRUE,
						xaxis = list(showticklabels = FALSE),
						yaxis = list(title = "RFU", range = c(0,20000))
					)  %>%
					config(edits = list(shapePosition = TRUE))
				
				event_register(p, "plotly_relayout")
				p
			})
			## 输出 r-square
			output[[paste0("model", i-1)]] <- renderText({
				usePeaks <- unname(attr(fsa, "ladderExact"))
				useSizes <- as.integer(names(attr(fsa, "ladderExact")))
				if (any(is.na(usePeaks))) {
					rSquared <- 0
				} else {
					model <- lm(useSizes ~ usePeaks)
					rSquared <- summary(model)$adj.r.squared
				}
				# 设置字体大小和边框样式
				if (rSquared < 0.999) {
					style <- "font-size: larger; color: red;"
					icon <- "&#9888;"
				} else {
					style <- "font-size: larger; color: green;"
					icon <- "&#10004;"
				}
				# 生成HTML代码
				HTML(paste("<div style='", style, "'>MQ: ", rSquared, " ", icon, "</div>"))
			})
		}
	})
}
## 重比对更新
fsa_updata <- function(fsa_copy, fsa_data, selected_colors, new_info, anno_df, input, output, session) {
	lapply(seq_len(length(fsa_copy) * 2), function(k) {
		if (k %% 2 == 0) {
			observeEvent(event_data("plotly_relayout", source = paste0("align", k)), {
				fsa <- fsa_copy[[k/2]]
				relayout_data <- event_data("plotly_relayout", source = paste0("align", k))
				# 获取所有影响到的 shape 的索引和 x0 值
				mod_indice <- grep("^shapes\\[[0-9]+\\]\\.x0$", names(relayout_data), value = TRUE)
				if (length(mod_indice) != 0) {
					mod_value <- unname(unlist(relayout_data[mod_indice]))
					for (j in seq_along(mod_indice)) {
						shape_index <- as.numeric(gsub("^shapes\\[([0-9]+)\\]\\.x0$", "\\1", mod_indice[j]))
						# 更新 ladderExact 值
						attr(fsa, "ladderExact")[shape_index + 1] <- mod_value[j]
					}

					usePeaks <- unname(attr(fsa, "ladderExact"))
					useSizes <- as.integer(names(attr(fsa, "ladderExact")))
					model <- lm(useSizes ~ usePeaks)
					rSquared <- summary(model)$adj.r.squared
					# Check alignment
					output[[paste0("model", k-1)]] <- renderText({
						# 设置字体大小和边框样式
						if (rSquared < 0.999) {
							style <- "font-size: larger; color: red;"
							icon <- "&#9888;"
						} else {
							style <- "font-size: larger; color: green;"
							icon <- "&#10004;"
						}
						# 生成HTML代码
						HTML(paste0("<div style='", style, "'>MQ: ", rSquared, " ", icon, "</div>"))
					})
					attr(fsa, "ladderModel") <- coefficients(model)
					attr(fsa,"ladderR2") <- rSquared
					fsa_copy[[k/2]] <<- fsa
				}
			})
			observeEvent(input$applyalign, {
				save_info(fsa_copy)
				ladder_plot(fsa_copy, new_info, anno_df, input, output, session)
				fsa_plot(fsa_copy, selected_colors, new_info, anno_df, input, output, session)
			})
			observeEvent(input$resetalign, {
				fsa_copy <<- fsa_data
				## 不重置到最初状态，而是重置这一次的realign操作
				ladder_plot(fsa_data, new_info, anno_df, input, output, session)
			})
		}
	})
}
## 调试rSquared
rsqu_out <- function(fsa_copy, info) {
	usePeaks <- unname(attr(fsa_copy[[1]], "ladderExact"))
	useSizes <- as.integer(names(attr(fsa_copy[[1]], "ladderExact")))
	model <- lm(useSizes ~ usePeaks)
	rSquared <- summary(model)$adj.r.squared
	print(info)
	print(rSquared)
}
## 重比对更新
save_info <- function(fsa_data) {
	peak_data <- data.frame()
	mapping_list = list()
	lapply(seq_along(fsa_data), function(i) {
		fsa <- fsa_data[[i]]
		temp_data <- data.frame()
		xcor <- attr(fsa, "ladderModel")[2] * 1:nrow(fsa) + attr(fsa, "ladderModel")[1]
		for (h in colnames(fsa)) {
			temp <- data.frame(sample = attr(fsa, "metaData")$sample,
								channel = h,
								xcor = xcor,
								fsa_value = fsa[, h])
			## 筛选peaks
			temp <- temp[temp$xcor > 30 & temp$fsa_value > 0, ]
			temp_data <- rbind(temp_data, temp)
		}
		peak_data <<- rbind(peak_data, temp_data)
		rSquared <- attr(fsa, "ladderR2")
		if (any(is.na(attr(fsa, "ladderExact")))) {
			mapping_info = c(attr(fsa, "metaData")$sample,0)
		} else {
			mapping_info = c(attr(fsa, "metaData")$sample,round(rSquared, 6))
		}
		mapping_list <<- append(mapping_list, mapping_info)
	})
	fwrite(peak_data, file = "www/peak.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
	mapping_df <- data.frame(matrix(unlist(mapping_list), ncol = 2, byrow = TRUE))
	colnames(mapping_df) <- c("Sample", "MQ")
	write.csv(mapping_df, file = "www/mapping_info.csv", row.names = FALSE)
}

library(DT)
library(shiny)
library(plotly)
library(shinyjs)
library(data.table)
library(shinyalert)
library(shinydashboard)
library(shinycssloaders)

fsa_data <- data.frame()

ui <- dashboardPage(
	title = 'ACMSI',
	dashboardHeader(
		title = tags$h3(
			style = "font-family: Archivo Black; font-size: 20px; color: white; margin-top: 12px; font-style: italic; font-weight: bold; text-align: center",
			"ACMSI"
		),
		tags$li(class = "dropdown",
				class = "dropdown-toggle",
				a(icon("dna"), "CE Fragment Analysis", href = "#", style = "color: #4D4D4D;")
		),
		dropdownMenuOutput("notification_menu")
	),
	dashboardSidebar(
		useShinyjs(),
		br(),
# 第一部分 - 标准品
		actionButton("side_std", label = HTML('<i class="fa-solid fa-flask"></i>&nbsp;&nbsp;&nbsp; Work Flow'), 
            style = "color: white; background-color: #444444; width: 200px; height: 40px; font-size: 1.2em; text-align: left;"),
		div(id = "sd_std",
			hidden = TRUE,
			uiOutput("std"),
			actionButton("view_std", label = HTML('<b>Edit</b> <i class="fa-solid fa-pen"></i>')),
		),
		br(),
# 第二部分 - 绘图
		actionButton("side_plot", label = HTML('<i class="fa-solid fa-chart-line"></i>&nbsp;&nbsp;&nbsp; Plot'), 
            style = "color: white; background-color: #444444; width: 200px; height: 40px; font-size: 1.2em; text-align: left;"),
		div(id = "sd_plot",
			hidden = TRUE,
			fileInput("annofile", "select annotation file (optional) ----- ## anno line ##", accept = c('.txt','.csv')),
			fileInput("infofile", "select sample info file (optional) ----- ## sort ##", accept = c('.txt','.csv')),
			fileInput("files", "select FSA files", multiple = TRUE, accept = c('.fsa')),
			actionButton("submit", label = HTML('<i class="fa-solid fa-sync fa-spin"></i> <b>Start!</b>')),
		),
		br(),
# 第三部分 - 辅助判断
		actionButton("side_classify", label = HTML('<i class="fa-solid fa-check-square"></i>&nbsp;&nbsp;&nbsp; Classification'), 
            style = "color: white; background-color: #444444; width: 200px; height: 40px; font-size: 1.2em; text-align: left;"),
		div(id = "sd_classify",
			hidden = TRUE,
			fileInput("anno", "select annotation file", multiple = FALSE, accept = c('.txt','.csv')),
			fileInput("info", "select sample info file", multiple = FALSE, accept = c('.txt','.csv')),
			actionButton("classify", label = HTML('<i class="fa-solid fa-sync fa-spin"></i> <b>GO!</b>'))
		)),
# 绘图区域
	dashboardBody(
		useShinyjs(),
		tabsetPanel(
			tabPanel("Plot", 
					br(),
					fluidRow(
					column(4, sliderInput("x_range", "X axis", min = -10, max = 300, value = c(-10, 300))),
					column(4, sliderInput("y_range", "Y axis", min = 0, max = 50000, value = c(0, 50000))),
					## 选择通道
					column(4, tags$style("
									.color-button {
									width: 70px;
									height: 30px;
									margin: 5px;
									border: 2px solid #aaa;
									border-radius: 5px;
									cursor: pointer;
									}
									.selected-button {
									border-color: #000000 !important;
									}
									.channel-label {
										font-size: 1em;
										font-weight: bold;
									}
								"),
								div(
									id = "color-buttons",
									h4("select channel", class = "channel-label"),
									actionButton("blue", "",
												class = "color-button",
												style = "background-color: #0066FF; border-color: #0066FF"),
									actionButton("green", "",
												class = "color-button",
												style = "background-color: #00FF00; border-color: #00FF00"),
									actionButton("yellow", "",
												class = "color-button",
												style = "background-color: #FFFF00; border-color: #FFFF00")
								)
							)
													),
					actionButton("view_mapping", label = HTML('<i class="fa-solid fa-list"></i> MQ Info'), 
								style = "color: #3c8dbc; background-color: white; width: 150px; height: 40px; font-size: 1.2em; text-align: middle; display:none"),
					shinycssloaders::withSpinner(uiOutput("plots"))),
			tabPanel("Calibration",
					br(),
					shinyjs::hidden(actionButton("applyalign", label = HTML('<i class="fa-solid fa-check"></i> apply'), 
								style = "color: #3c8dbc; background-color: white; width: 150px; height: 40px; font-size: 1.2em; text-align: middle")),
					shinyjs::hidden(actionButton("resetalign", label = HTML('<i class="fa-solid fa-undo"></i> reset'), 
								style = "color: #3c8dbc; background-color: white; width: 150px; height: 40px; font-size: 1.2em; text-align: middle")),
					shinycssloaders::withSpinner(uiOutput("realigns"))),
			tabPanel("Classification", 
					br(),
					actionButton("view_quality", label = HTML('<i class="fa-solid fa-list"></i> QC Info'), 
								style = "color: #3c8dbc; background-color: white; width: 150px; height: 40px; font-size: 1.2em; text-align: middle; display:none"),
					br(),
					br(),
					shinyjs::hidden(htmlOutput("classify_result")),
					shinycssloaders::withSpinner(imageOutput("image"))),
		),
		useShinyjs()
	)
)

server <- function(input, output, session) {
## 读取fsa
	hideSpinner("image")
	fsa_data_reactive <- reactiveVal()
	## 通道选择
	selected_colors <- reactiveVal(c("blue", "green", "yellow"))
	shinyjs::toggleClass("green", "selected-button")
	shinyjs::toggleClass("blue", "selected-button")
	shinyjs::toggleClass("yellow", "selected-button")
	## 弹窗信息
	notification_data <- reactiveValues(
		message_type = NULL,
		message_content = NULL,
		message_icon = NULL
	)
	## 对应dye和颜色的Database
	anno_df <- data.frame(
		Channel = c("NED", "VIC", "6-FAM", "PET", "dR110", "dR6G", "dTAMRA", "LIZ", "dROX", "5-FAM", "JOE", "ROX", "HEX", "TAMRA"),
		Color = c("yellow", "green", "blue", "red", "blue", "green", "yellow", "orange", "red", "blue", "green", "red", "green", "red")
	)

	observeEvent(input$submit, {
		if (is.null(input$std) || nchar(input$std) == 0) {
			shinyalert("Please select STD first！", type = "warning")
		}
		else if (is.null(input$files)) {
			shinyalert("Please select FSA files！", type = "warning")
		}
		else {
			showSpinner("plots")
			fsa_data <- fsa_convert(input$std, input$files, anno_df, session)
			if ("stdNotFound" %in% fsa_data) {
				fsa <- read.fsa(input$files$datapath[[1]])
				hideSpinner("plots")
				shinyjs::enable("submit")
				shinyalert(paste("STD ",input$std, " channels are not available ( Suggestions: ", paste(colnames(fsa), collapse=", "), ")"), type = "error")
			} else {
			fsa_data_reactive(fsa_data)
			existing_info <- read.table("www/mapping_info.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
			if (any(existing_info$MQ == 0)) {
				shinyalert(paste("Check ",input$std," channels again! Can't find CE bands."), type = "error")
			}
			danger_rows <- existing_info[(existing_info[,2] < 0.999), 'Sample']
			if (length(danger_rows) == 0) {
				message_content <- "Align: all samples align successfullt"
				notification_data$message_type <- "success"
				notification_data$message_icon <- shiny::icon("check")
			} else {
				message_content <- paste("Align： the following samples failed in aligning (MQ < 0.999), not suggested for auto-classify: ", paste(danger_rows, collapse = ", "))
				notification_data$message_type <- "danger"
				notification_data$message_icon <- shiny::icon("exclamation-triangle")
			}
			notification_data$message_content <- paste(notification_data$message_content, message_content, sep = "<br>")
			shinyjs::show("view_mapping")
			showSpinner("plots")
			}
		}
	})
## 绘制
    observe({
		fsa_data <- fsa_data_reactive()
		selected_colors <- selected_colors()
		new_info <- data.frame(name = character(), Sample = character(), stringsAsFactors = FALSE)
		for (i in 1:length(fsa_data)) {
		new_info <- rbind(new_info, data.frame(name = "", Sample = ""))
		}
        if (!is.null(fsa_data)) {
			if (!is.null(input$infofile)) {
				info_df <- read.table(input$infofile$datapath, header = TRUE, sep = ",", stringsAsFactors = FALSE)
				new_info <- data.frame(name = character(), Sample = character(), stringsAsFactors = FALSE)
				for (i in 1:nrow(info_df)) {
					new_info <- rbind(new_info, data.frame(name = c(info_df$Tname[i], info_df$Nname[i]), Sample = info_df$Sample[i]))
				}
				new_info$name <- as.character(new_info$name)
				sorted_fsa_data <- fsa_data[match(new_info$name, sapply(fsa_data, function(x) attr(x, "metaData")$sample))]
				fsa_plot(sorted_fsa_data, selected_colors, new_info, anno_df, input, output, session)
				ladder_plot(sorted_fsa_data, new_info, anno_df, input, output, session)
				fsa_updata(sorted_fsa_data, sorted_fsa_data, selected_colors, new_info, anno_df, input, output, session)
			} else {
				fsa_plot(fsa_data, selected_colors, new_info, anno_df, input, output, session)
				ladder_plot(fsa_data, new_info, anno_df, input, output, session)
				fsa_updata(fsa_data, fsa_data, selected_colors, new_info, anno_df, input, output, session)
			}
			## 重比对更新状态
			
			shinyjs::show("applyalign")
			shinyjs::show("resetalign")
        }
    })

	## 主要绘图区
    output$plots <- renderUI({
		fsa_data <- fsa_data_reactive()
        if (!is.null(fsa_data)) {
			plot_output_list <- lapply(seq_along(fsa_data), function(i) {
				plotlyOutput(outputId = paste0("plot", i))
			})
			do.call(tagList, plot_output_list)
        }
    })

	## 重比对
	output$realigns <- renderUI({
		fsa_data <- fsa_data_reactive()
		if (!is.null(fsa_data)) {
			realign_output_list <- lapply(seq_len(length(fsa_data) * 2), function(i) {
			if (i %% 2 == 1) {  # 输出文字框
				htmlOutput(outputId = paste0("model", i))
			} else {  			# 输出绘图框
				plotlyOutput(outputId = paste0("realign", i))
			}
			})
			do.call(tagList, realign_output_list)
		}
	})

	## 自动判别
    observeEvent(input$classify, {
		if (is.null(input$anno)) {
			shinyalert("Please select annotation file of panels first!", type = "warning")
		} 
		else if (is.null(input$info)) {
			shinyalert("Please selection the information files of samples first!", type = "warning")
		}
		else {
			shinyjs::disable("classify")
			showSpinner(id="image")

			peak_path <- file.path(getwd(), "www", "peak.txt")
			command <- paste(".\\dist\\MSIclassify\\MSIclassify.exe -i", input$info$datapath, " -p ", peak_path, " -a ", input$anno$datapath)
			cat(command, "\n")
			# 执行命令
			system(command, wait = TRUE)

			output$image <- renderImage({
				list(src = "www/plot_all_img.png")
			}, deleteFile = FALSE)
			shinyjs::enable("classify")
			shinyjs::show("classify_result")
			shinyjs::show("view_quality")	
		}
## 小弹窗提示
		existing_info <- read.table("www/quality_info.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
		colnames(existing_info) <- c("Sample", "Marker", "MaxRFU", "shiftBP")
		existing_info$QC <- ifelse(existing_info[, 3] < 1000, "Failed",
									ifelse(existing_info[, 3] < 3000, "Warning", "PASS"))

		warning_rows <- existing_info[existing_info$QC %in% c("Warning"), c("Sample", "Marker")]
		danger_rows <- existing_info[existing_info$QC %in% c("Failed"), c("Sample", "Marker")]
		if ((length(warning_rows$Sample) > 0) && (length(danger_rows$Sample) == 0)) {
		message_content <- paste("<br>QC： the followings' qc are at risk (RFU < 3000): ", paste(paste(warning_rows$Sample, warning_rows$Marker, sep = "-"),  sep = "<br>"), sep = "<br>")
		notification_data$message_type <- "warning"
		notification_data$message_icon <- shiny::icon("warning")
		}
		else if ((length(warning_rows$Sample) == 0) && (length(danger_rows$Sample) > 0)) {
			message_content <- paste("<br>QC： the followings' qc failed (RFU < 1000): ", paste(paste(danger_rows$Sample, danger_rows$Marker, sep = "-"), sep = "<br>"), sep = "<br>")
			notification_data$message_type <- "danger"
			notification_data$message_icon <- shiny::icon("exclamation-triangle")
		} 
		else if ((length(warning_rows$Sample) > 0) && (length(danger_rows$Sample) > 0)) {
			message_content <- paste("<br>QC： the followings' qc failed (RFU < 1000): ", paste(paste(danger_rows$Sample, danger_rows$Marker, sep = "-"), sep = "<br>"), "<br> the followings' qc is are at risk  (RFU < 3000): ", paste(paste(warning_rows$Sample, warning_rows$Marker, sep = "-"), sep = "<br>"), sep = "<br>")
			notification_data$message_type <- "danger"
			notification_data$message_icon <- shiny::icon("exclamation-triangle")
		}
		else {
			message_content <- "<br>QC: all markers' qc passed"
			notification_data$message_type <- "success"
			notification_data$message_icon <- shiny::icon("check")
		}
		notification_data$message_content <- paste(notification_data$message_content, message_content, sep = "<br>")
		})
		
	output$notification_menu <- renderMenu({
		dropdownMenu(type = "notification", badgeStatus = notification_data$message_type,
				menuItem(div(style = "overflow: auto;", HTML(notification_data$message_content)), icon = notification_data$message_icon))
	})


## 质控信息
	observeEvent(input$view_quality, {
		if (file.exists("www/quality_info.csv") && file.size("www/quality_info.csv") > 0) {
			output$quality_info_table <- DT::renderDataTable({
				existing_info <- read.table("www/quality_info.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
				existing_info <- existing_info[, -c(ncol(existing_info)-1, ncol(existing_info))]
				colnames(existing_info) <- c("Sample", "Marker", "MaxRFU","shiftBP")
				existing_info$QC <- ifelse(existing_info[, 3] < 1000, "Failed",
											ifelse(existing_info[, 3] < 3000, "Warning", "PASS"))
				datatable(existing_info, 
							options = list(
							dom = 'Bfrtip',
							paging = TRUE
							)) %>%
					formatStyle(
					columns = "QC",
					valueColumns = "QC",
					target = "row",
					backgroundColor = styleEqual(c("Failed", "Warning"), c("red", "yellow"))
				)
			})
			showModal(modalDialog(
			title = "QC",
			div(style = "max-height: 500px; overflow: auto;", 
				DT::dataTableOutput("quality_info_table")),
			footer = tagList(
				modalButton("Close")
			)
			))
		} else {
			shinyalert("Please run automated classification first", type = "warning")
		}
	})
## 自动判别
	output$classify_result <- renderText({
		classify_data <- read.csv("www/quality_info.csv")
		# 筛选mark等于1的行
		filtered_data <- subset(classify_data, mark == 1)
		output_text <- paste("<i class='fas fa-exclamation-triangle' style='color: red;'></i> <b>Sample ", filtered_data$sample, ": Marker", filtered_data$site, "possibly be MSI</b>", collapse = "<br/>")
		output_text <- paste("<span style='color: red;'>", output_text, "</span>")
		
		HTML(output_text)
	})
## 比对信息
	observeEvent(input$view_mapping, {
		if (file.exists("www/mapping_info.csv") && file.size("www/mapping_info.csv") > 0) {
			## 点击后生成一个图表框，保持实时更新
			output$mapping_info_table <- DT::renderDataTable({
				existing_info <- read.table("www/mapping_info.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
				colnames(existing_info) <- c("Sample", "MQ(>0.999)")
				existing_info$MQ <- ifelse(existing_info[, 2] < 0.999, "Failed", "PASS")
				datatable(existing_info, 
							options = list(
							dom = 'Bfrtip',
							paging = TRUE
							)) %>%
					formatStyle(
					columns = "MQ",
					valueColumns = "MQ",
					target = "row",
					backgroundColor = styleEqual(c("Failed"), c("red"))
				)
			})
			showModal(modalDialog(
			title = "size calling",
			div(style = "max-height: 500px; overflow: auto;", 
				DT::dataTableOutput("mapping_info_table")),
			footer = tagList(
				modalButton("Close")
			)
			))
		} else {
			shinyalert("Please run plot first!", type = "warning")
		}
	})

## 选择标准品
	output$std <- renderUI({
		if (file.exists("www/std_info.txt") && file.size("www/std_info.txt") > 0) {
			existing_info <- read.table("www/std_info.txt", header = TRUE, sep = "\t")
			selectInput("std", "select STD", choices = existing_info$name, selected = NULL)
		} else {
			shinyalert("No saved STD infomation, please add by `Edit` first!", type = "warning")
			selectInput("std", "select STD", choices = NULL, selected = NULL)
		}
	})
## 查看标准品明细
	observeEvent(input$view_std, {
		if (file.exists("www/std_info.txt") && file.size("www/std_info.txt") > 0) {
			output$std_info_table <- DT::renderDataTable({
				existing_info <- read.table("www/std_info.txt", header = TRUE, stringsAsFactors = FALSE)
				colnames(existing_info) <- c("STD", "Channel", "Size")
				datatable(existing_info, 
					editable = TRUE,
					selection = 'multiple',
					options = list(
					dom = 'Bfrtip',
					paging = TRUE
					)
				)
			})
			showModal(modalDialog(
			title = "STD infomation",
			div(style = "max-height: 500px; overflow: auto;", 
				DT::dataTableOutput("std_info_table")),
			actionButton("add_btn", "ADD"),
			actionButton("delete_btn", "Delete"),
			footer = tagList(
				modalButton("Close")
			)
			))
		} else {
			shinyalert("No saved STD infomation, please add by `Edit` first!", type = "warning")
		}
		})

	observeEvent(input$delete_btn, {
		selected_rows <- as.numeric(input$std_info_table_rows_selected)
		if(length(selected_rows) > 0){
			showModal(modalDialog(
			title = "Delete Still",
			"Still want to delete the selected STD ?",
			footer = tagList(
				modalButton("Cancel"),
				actionButton("confirm_delete", "Delete Still", class = "btn-danger")
			)
			))
		} else {
			shinyalert("Please select the STD to be deleted!", type = "warning")
		}
	})
	observeEvent(input$confirm_delete, {
		selected_rows <- as.numeric(input$std_info_table_rows_selected)
		existing_info <- read.table("www/std_info.txt", header = TRUE, stringsAsFactors = FALSE)
		existing_info <- existing_info[-selected_rows,]
		write.table(existing_info, "www/std_info.txt", sep = "\t", row.names = FALSE)
		# 更新select
		output$std <- renderUI({
			existing_info <- read.table("www/std_info.txt", header = TRUE, sep = "\t")
			selectInput("std", "select STD", choices = existing_info$name, selected = NULL)
		})
		removeModal()
		shinyalert(paste("STD", existing_info[selected_rows]$name, "has already been deleted"), type = "success")
	})
## 添加标准品
	observeEvent(input$add_btn, {
		showModal(modalDialog(
			title = "ADD STD",
			textInput("name", "STD name"),
			textInput("channel", "STD channel"),
			textInput("ladder", "STD size range"),
			actionButton("stdsubmit", "submit"),
			footer = tagList(
			modalButton("cancel")
			)
		))
	})

	observeEvent(input$stdsubmit, {
			if (!is.null(input$name) && (nchar(input$name) != 0) &&
				!is.null(input$channel) && (nchar(input$channel) != 0) &&
				!is.null(input$ladder)  && (nchar(input$ladder) != 0)) {
				if (!grepl("^(\\d{1,4},)+\\d{1,4}$", input$ladder)) { 
					shinyalert("Ladder must be a comma-separated combination of numbers, such as 20,40,60. Check the format and enter a ladder again", type = "error")
				} 
				else {
					if (file.exists("www/std_info.txt")) {
						existing_info <- read.table("www/std_info.txt", header = TRUE, sep = "\t")
						if (input$name %in% existing_info$name) {
							shinyalert("STD name already exists, please re-enter!", type = "warning")
						} else {
							std_info <- data.frame(name = input$name, channel = input$channel, ladder = input$ladder)
							new_info <- rbind(existing_info, std_info)
							write.table(new_info, "www/std_info.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)
							updateSelectInput(session, "std", choices = new_info$name)
							shinyjs::hide("new_inputs")
							shinyalert(paste("STD", std_info$name, "is added successfully!"), type = "success")
							output$std_info_table <- DT::renderDataTable({
								existing_info <- read.table("www/std_info.txt", header = TRUE, stringsAsFactors = FALSE)
								colnames(existing_info) <- c("Sample", "Channel", "Size")
								datatable(existing_info, 
											editable = TRUE,
											selection = 'multiple',
											options = list(
											dom = 'Bfrtip',
											paging = TRUE
											)
								)
							})
							removeModal()
						}
					} else {
						std_info <- data.frame(name = input$name, channel = input$channel, ladder = input$ladder)
						write.table(std_info, "www/std_info.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)
						updateSelectInput(session, "std", choices = std_info$name)
						shinyjs::hide("new_inputs")
						shinyalert(paste("STD", std_info$name, "is added successfully!"), type = "success")
						output$std_info_table <- DT::renderDataTable({
							existing_info <- read.table("www/std_info.txt", header = TRUE, stringsAsFactors = FALSE)
							colnames(existing_info) <- c("STD", "Channel", "Size")
							datatable(existing_info, 
										editable = TRUE,
										selection = 'multiple',
										options = list(
										dom = 'Bfrtip',
										paging = TRUE
										)
							)
						})
						removeModal()
					}
				}
			}
			else {
				shinyalert("Remain STD info are not provided, check again!", type = "error")
			}
	})
## 通道选择
	observeEvent(input$green, {
		if (input$green %% 2 == 0) {
		shinyjs::toggleClass("green", "selected-button")
		selected_colors(c(selected_colors(), "green"))
		} else {
		shinyjs::toggleClass("green", "selected-button")
		selected_colors(selected_colors()[selected_colors() != "green"])
		}
	})
	
	observeEvent(input$blue, {
		if (input$blue %% 2 == 0) {
		shinyjs::toggleClass("blue", "selected-button")
		selected_colors(c(selected_colors(), "blue"))
		} else {
		shinyjs::toggleClass("blue", "selected-button")
		selected_colors(selected_colors()[selected_colors() != "blue"])
		}
	})
	
	observeEvent(input$yellow, {
		if (input$yellow %% 2 == 0) {
		shinyjs::toggleClass("yellow", "selected-button")
		selected_colors(c(selected_colors(), "yellow"))
		} else {
		shinyjs::toggleClass("yellow", "selected-button")
		selected_colors(selected_colors()[selected_colors() != "yellow"])
		}
	})
## 按钮收放
	observeEvent(input$side_std, {
		shinyjs::toggle("sd_std")
	})
		observeEvent(input$side_plot, {
		shinyjs::toggle("sd_plot")
	})
		observeEvent(input$side_classify, {
		shinyjs::toggle("sd_classify")
	})
}

# 运行Shiny应用程序
shinyApp(ui, server)
# shinyApp(ui, server, options = list(host = "0.0.0.0", port = 4878))
