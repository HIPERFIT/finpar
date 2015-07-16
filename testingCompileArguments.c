// THIS FILE CAN BE COMPILED WITH: 
// gcc -o testingCompileArguments testingCompileArguments.c -D STRING='"Hello World"'
// This tells the preprocessor to replace STRING by the value on the RHS of the equality sign.
#include <stdio.h>

int main (int argc, char *argv[]) {
  printf("String is: %s\n", STRING);
  return 0;
}
