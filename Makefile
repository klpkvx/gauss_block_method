FLAGS = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
OPT = -O3
CC = g++

all: a.out

%.o: %.cpp
	${CC} -c ${FLAGS} ${OPT}  $<

%.out: %.o init.o matrix.o solution.o
	${CC} $^ -o  $@

clean:
	rm -f *.out *.o *.d