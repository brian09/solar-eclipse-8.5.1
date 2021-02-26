EXECS = ibdprep ibdmat ibdmc allfreq relate getmeans multipnt mrgibd makeped unknown genfreq #dommsibs mkibdf
all: $(EXECS)

# Mac
# Note: -finit-local-zero not available
FC ?= gfortran
CC ?= gcc
SAFELIB_PATH ?= /usr/local
CFLAGS = -static -O2 -I$(SAFELIB_PATH)/include
FFLAGS = -O2
LDFLAGS = -static  -O2
ifneq ($(OSTYPE), linux-gnu)
CFLAGS =  -O2 -I$(SAFELIB_PATH)/include
LDFLAGS = -O2
endif
NOVECTOR = 
SUNMATH =
INSTALLER = cp
INSTALLER2 = $(INSTALL_PATH)

ibdprep: ibdprep.o
	$(CC) $(LDFLAGS) -o ibdprep ibdprep.o -L$(SAFELIB_PATH)/lib -lsafe -lm

ibdmat: ibdmat.o
	$(CC) $(LDFLAGS) -o ibdmat ibdmat.o

ibdmc: ibdmc.o
	$(FC) $(LDFLAGS) -o ibdmc ibdmc.o

allfreq: allfreq.o fallfrq.f
	$(FC) $(LDFLAGS)  $(NOVECTOR) -o allfreq allfreq.o fallfrq.f

genfreq: allfreq.o fgenfrq.f
	$(FC) $(LDFLAGS) $(NOVECTOR) -o genfreq allfreq.o fgenfrq.f

relate: relate.o
	$(CC) $(LDFLAGS) -o relate relate.o

getmeans: getmeans.o
	$(CC) $(LDFLAGS) -o getmeans getmeans.o -L$(SAFELIB_PATH)/lib -lsafe -lm

multipnt: dgedifa.o multipnt.o
	$(CC) $(LDFLAGS) -o multipnt multipnt.o dgedifa.o -L$(SAFELIB_PATH)/lib -lsafe $(SUNMATH) -lm

mrgibd: mrgibd.o
	$(CC) $(LDFLAGS) -o mrgibd mrgibd.o -L$(SAFELIB_PATH)/lib -lsafe

mkibdf: mkibdf.o
	$(CC) $(LDFLAGS) -o mkibdf mkibdf.o -L$(SAFELIB_PATH)/lib -lsafe

makeped: makeped.c
	$(CC) -o makeped makeped.c

unknown: unknown.o
	$(CC) $(LDFLAGS) -o unknown unknown.o

install:
	for i in $(EXECS) ; do ( $(INSTALLER) $$i $(INSTALLER2) ) ; done

test:
	for i in $(EXECS) ; do ( /usr/sbin/install -f ../newbin -o $$i ) ; done

clean:
	rm -f $(EXECS) *.o
