# <In terms of Make, a phony target is simply a target that is 
# always out - of - date, so whenever you ask make < phony_target>,
# it will run, independent from the state of the file system>

.PHONY: clean
CC = g++
CFLAGS = -Wall -Wextra -Werror -std=c++17
COVERAGE_FLAGS = -fprofile-arcs -ftest-coverage
LIBS = -lgtest -lpthread
MAIN_FILE = s21_matrix_oop.cc
TESTS_FILE = unit_tests.cc

all: s21_matrix_oop.a

s21_matrix_oop.a: s21_matrix_oop.o
	ar rc s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

s21_matrix_oop.o:
	$(CC) $(CFLAGS) -c $(COVERAGE_FLAGS) $(MAIN_FILE) -o s21_matrix_oop.o $(LIBS)

test: s21_matrix_oop.o
	$(CC) $(CFLAGS) $(COVERAGE_FLAGS) s21_matrix_oop.o $(TESTS_FILE) -o test.out $(LIBS)
	./test.out

#--no-external flag is for not include any files, except our src
gcov_report: test
	lcov -t "./test" -o test.info --no-external -c -d ./
	genhtml -o report test.info
	open report/index.html
	rm *.info

leaks: test
	leaks -atExit -- ./test.out

style:
	cp ../materials/linters/.clang-format .
	clang-format -n *.cc *.h

style-fix:
	cp ../materials/linters/.clang-format .
	clang-format -i *.cc *.h

clean:
	rm -f *.a
	rm -f *.o
	rm -f *.out
	rm -f *.gcno
	rm -f *.gcda
	rm -rf test.out.dSYM
	rm -rf report