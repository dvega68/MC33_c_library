# Project: libMC33
NAME     = MC33
!IFNDEF MACHINE
MACHINE  = x64
!ENDIF
LIBDIR   = lib\$(MACHINE)
SRC      = source\libMC33.c
OBJ      = libMC33.obj
BIN      = $(LIBDIR)\$(NAME).lib
CFLAGS   = /GL /W3 /Gy /I".\include" /O2 /Zc:inline /fp:fast /Gd /Oi /MD /EHsc /nologo /Ot /D _CRT_SECURE_NO_WARNINGS
LIBFLAGS = /LTCG /MACHINE:$(MACHINE) /NOLOGO
 
all: $(BIN)

$(OBJ): $(SRC)
	$(CC) $(CFLAGS) /c $(SRC)

$(BIN): $(OBJ)
	@if not exist $(LIBDIR) mkdir $(LIBDIR)
	LIB /OUT:$(BIN) $(LIBFLAGS) $(OBJ)

clean:
	-1 del $(OBJ)
	-1 del $(BIN)
