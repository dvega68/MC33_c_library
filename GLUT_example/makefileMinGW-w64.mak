#--------------------------------------------------------
CC       = gcc
SOURCE   = TestMC33_glut.c
OBJ      = $(SOURCE:.c=.o)
LIBS     = -static -lfreeglut -lglu32 -lopengl32 -lgdi32 -lwinmm -lm -m64 -s
CINCS    =
BIN      = TestMC33_glut
OPTIM		 = -Ofast -m64 -Wall -Wextra -funroll-loops
CFLAGS   = -Wno-unused-parameter $(CINCS) -DFREEGLUT_STATIC $(OPTIM)
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

all:	all-before $(BIN) all-after

$(BIN): $(OBJ)
	#$(SOURCE)
	$(CC) -o $(BIN) $(OBJ) $(LIBS)

clean: clean-custom
	$(RM) $(OBJ) $(BIN)
#--------------------------------------------------------
