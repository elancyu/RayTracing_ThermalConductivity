LIBS = -lm
CFLAGS = -Wall -O3 -g -I/usr/local/include
LFLAGS = -L/HOME/jyuaq/RayTracingC/RayTracing
LD = cc $(LFLAGS)
CC = gcc -c

OBJS = main.o\
       BuildGeometry.o\
       Calculation.o\
	   FreeMem.o\
	   GlobalVariables.o\
	   Rand.c\
	   PerformRayTracing.o\
	   Transmission.o
	   
	   
main: $(OBJS)
	$(LD) -o RayTracing $(OBJS) $(LIBS)
	 
.c.o:
	$(CC) $(CFLAGS) $<

>PHONY: clean
clean:
	rm ${OBJS} *.dat *.log