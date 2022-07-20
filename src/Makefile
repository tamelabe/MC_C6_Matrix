OS := $(shell uname -s)
CC = gcc
FLAG = -Wall -Werror -Wextra
FLAG_C = -c -Wall -Werror -Wextra
CL_11 = -std=c11
ALL_FILE = s21_matrix.c
ALL_FILE_O = s21_matrix.o

ifeq ($(OS), Darwin)
	FLAGS_PLATFORM = -lcheck
else
	FLAGS_PLATFORM = -lcheck -lsubunit -lrt -lm -lpthread
endif

all: clean s21_matrix.a

s21_matrix.a:
	@$(CC) $(FLAG_C) $(ALL_FILE)
	@ar rcs s21_matrix.a $(ALL_FILE_O)

test: s21_matrix.a
	@gcc tests.c s21_matrix.a $(FLAGS_PLATFORM) -g -o tests.o
	@./tests.o

gcov_report:
	@$(CC) --coverage $(ALL_FILE) tests.c $(FLAGS_PLATFORM) -o gcov_report.o
	@./gcov_report.o
	@lcov -t s21_matrix_tests -o s21_matrix_tests.info -c -d .
	@genhtml -o report s21_matrix_tests.info
	@open report/index.html

clean:
	@rm -f *.a
	@rm -f s21_matrix
	@rm -f *.o
	@rm -f *.txt
	@rm -f *.gcda
	@rm -f *.gcdo
	@rm -f *.gcno
	@rm -f *.info
	@rm -rf ./report

leaks: $(ALL_FILE_O)
	@CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./tests.o

style: $(ALL_FILE)
	@cp ../materials/linters/CPPLINT.cfg ./
	@python3 ./../materials/linters/cpplint.py --extensions=c $(ALL_FILE) tests.c
	@rm -f CPPLINT.cfg


check: style cppcheck