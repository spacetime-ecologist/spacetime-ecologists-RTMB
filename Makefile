#
#       run and save all opt$objective
# from each '_RTMB.R' file in this repository
#           and store in RTMB.res 
#

# define the directory as the current directory
RTMBDIR := $(CURDIR)

# find _RTMB.R files in the current directory, ignore ch 10 
infiles := $(shell find $(RTMBDIR) -type f -name \*_RTMB.R \)

outfiles := $(shell echo $(infiles) | sed s/_RTMB.R/_RTMB.chk/g)

.PHONY: all clean list res

# main target
all: res

clean:
	@rm -f $(outfiles) RTMB.res

# write each opt$objective to a .chk file
$(CURDIR)/%.chk: $(CURDIR)/%.R
	echo 'source("$<", chdir=TRUE);' \
		'opts <- ls(pattern="^opt");' \
		'write(sapply(opts, function(x) paste(get(x)$$objective, collapse="\n")), file="$@");' | R --vanilla --slave

list:
	@echo $(infiles)
	@echo $(outfiles)

# combine .chk files into RTMB.res, delete indv .chk files
res: $(outfiles)
	@for f in $(outfiles); do \
		fname=$$(basename $$f .chk).R; \
		echo "$$fname" >> RTMB.res; \
		cat $$f >> RTMB.res; \
		echo "" >> RTMB.res; \
		rm $$f; \
	done
	@echo "Results collected in RTMB.res"
