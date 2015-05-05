FC	= ifort
FFLAGS	= -cm -w -02
LIBS	= 

COMPILE = $(FC) $(FFLAGS)
LINK	= $(FC) $(LDFLAGS)

TARGET	= gradmin
OBJS	= calcgrad.o chicalc.o f90lib.o Gradiant.o writeinput.o

$(TARGET): $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)
	
%.o:%.f90
	$(COMPILE) -c $<