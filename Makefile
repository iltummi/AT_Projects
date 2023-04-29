#---------------------------------------------------
# Target file to be compiled by default
#---------------------------------------------------
MAIN	=	first_main
#---------------------------------------------------
# CC is the compiler to be used
#---------------------------------------------------
CC	=	gcc
#---------------------------------------------------
# CFLAGS are the options passed to the compiler
#---------------------------------------------------
CFLAGS	=	-Wall
#---------------------------------------------------
# OBJS are the object files to be linked
#---------------------------------------------------
OBJ1	=	task
OBJ2	= 	predittore
OBJ3	=	motor
OBJ4	=	time_operations
OBJ5	=	camera_unico
OBJS	=	$(MAIN).o $(OBJ1).o $(OBJ2).o $(OBJ3).o $(OBJ4).o $(OBJ5).o
#---------------------------------------------------
# LIBS are the external libraries to be used
#---------------------------------------------------
LIBS	=	-lpthread -lwiringPi `allegro-config --libs` #-lsemaphore –lrt -lm
#---------------------------------------------------
# Dependencies
#---------------------------------------------------
$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

$(MAIN).o: $(MAIN).c
	$(CC) $(CFLAGS) -c $(MAIN).c

$(OBJ1).o: $(OBJ1).c
	$(CC) $(CFLAGS) -c $(OBJ1).c

$(OBJ2).o: $(OBJ2).c
	$(CC) $(CFLAGS) -c $(OBJ2).c

$(OBJ3).o: $(OBJ3).c
	$(CC) $(CFLAGS) -c $(OBJ3).c

$(OBJ4).o: $(OBJ4).c
	$(CC) $(CFLAGS) -c $(OBJ4).c

$(OBJ5).o: $(OBJ5).c
	$(CC) $(CFLAGS) -c $(OBJ5).c
#---------------------------------------------------
# Command that can be specified inline: make clean
#---------------------------------------------------
clean:
	rm -rf *o $(MAIN)