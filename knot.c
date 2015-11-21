// C program to input an ascii picture of a tangle
// and output the matrix it would compute in a Fibonacci quantum computer.
//
// We work over the field of integers mod 521, with q=5, so q^5 = -1.
// Basis vectors are enumerated by binary numbers with no consecutive 1's.
// 
// Input:
// |, \, / mean a vertical strand
// % means a positive crossing
// S means a negative crossing
// u means a cup
// n means a cap
// . means end
// spaces are ignored
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MOD 521
#define Q 5

#define QQ ((Q * Q) % MOD)
#define QQQ ((Q * QQ) % MOD)
#define QQQQ ((Q * QQQ) % MOD)
#define PHI ((Q + MOD - QQQQ) % MOD)
#define PHI_INV (PHI - 1)

#define CUP {1,0,  0,0,  PHI_INV,0,  0,0,  0,1}
#define CAP {1,0,1,0,0,  0,0,0,0,PHI}
#define UNCROSS {MOD-PHI_INV,0,QQQ,0,0,    0,MOD-QQ,0,0,0,    (QQQ*PHI_INV)%MOD,0,((MOD-Q)*PHI_INV)%MOD,0,0,    0,0,0,MOD-QQ,0,    0,0,0,0,QQQQ}
#define MERGE {MOD-PHI_INV,0,1,0,0,    0,1,0,0,0,    0,0,0,1,0}


#define BIFMAX 20 //largest fibonacci number we'll ever use
#define MAX 6765  //MAX = fib(BIFMAX)
#define MAXSTRING 1000 // longest string input

int fib[BIFMAX + 1];  // fibonacci numbers
int fibword[MAX];
  // the infinite fibonacci word (A003849 in OEIS).
  // The basis vectors are indexed by Fibbinary numbers
  // (ie. strings of bits that have no adjacent ones, see A003714)
  // and fibword[i] is the last bit of the ith basis vector.


int * mat[256];
int bif_rows[256];
int bif_cols[256];

void freeall()
{
  int i;
  for(i = 0; i < 256; i++)
  {
    if(mat[i] != NULL)
    {
      free(mat[i]);
    }
  }
}

void die(char* s)
{
  freeall();
  printf("ERROR: %s\n",s);
  exit(0);
}

int * copy_array(int x[], int size)
{
  int i;
  int * y;
  y = malloc(size * sizeof(int));
  for(i = 0; i < size; i++)
  {
    y[i] = x[i];
  }
  return(y);
}

void initialize()
{
  int i;
  int j;

  int cap[10] = 
          {1,0,1,0,0,
           0,0,0,0,PHI};
  int cup[10] =
          {1,        0,
           0,        0,
           PHI_INV,  0,
           0,        0,
           0,        1};
  int cross[25] = 
          {MOD-PHI_INV,            0,    MOD-QQ,             0,   0,
           0,                      QQQ,  0,                  0,   0,
           ((MOD-QQ)*PHI_INV)%MOD, 0,    (QQQQ*PHI_INV)%MOD, 0,   0,
           0,                      0,    0,                  QQQ, 0,
           0,                      0,    0,                  0,   MOD-Q};
  int id[4] = 
          {1,0,
           0,1};
  int split[15] =
          {MOD-PHI_INV,  0,  0,
           0,            1,  0,
           PHI_INV,      0,  0,
           0,            0,  1,
           0,            0,  0};

  //initialize mat[]
  for(i = 0; i < 256; i++)
  {
    mat[i] = NULL;
  }
  mat['^'] = copy_array(cap,10);
  bif_rows['^'] = 3;
  bif_cols['^'] = 5;
  mat['_'] = copy_array(cup,10);
  bif_rows['_'] = 5;
  bif_cols['_'] = 3;
  mat['%'] = copy_array(cross,25);
  bif_rows['%'] = 5;
  bif_cols['%'] = 5;
  mat['|'] = copy_array(id,4);
  bif_rows['|'] = 4;
  bif_cols['|'] = 4;
  mat['4'] = copy_array(split,20);
  bif_rows['4'] = 4;
  bif_cols['4'] = 5;

  //initialize fib[]
  fib[0] = 0;
  fib[1] = 1;
  for(i = 2; i <= BIFMAX; i++)
  {
    fib[i] = fib[i-1] + fib[i-2];
  }

  // Initialize fibword to the "fibonacci word" - see the comment when it was declared.
  fibword[0] = 0;
  fibword[1] = 1;
  for(i = 2; i < BIFMAX - 1; i++)
  {
    // Make a copy of the first fib[i] bits shifted by fib[i+1]
    for(j = 0; j < fib[i]; j++)
    {
      if(fib[i + 1] + j >= MAX)
      {
        die("error initializing fibword");
      }
      fibword[fib[i + 1] + j] = fibword[j];
    }
  }
}

void prettyprint()
{
  int i = 0;
  int j = 0;
  int * A = mat['!'];
  int r = fib[bif_rows['!']];
  int c = fib[bif_cols['!']];

  if(A == NULL)
  {
    printf("Tried to print a null pointer.\n");
    return;
  }
  if(r * c > 1000)
  {
    printf("%d by %d is too big to pretty print.\n", r, c);
    return;
  }

  for(i = 0; i < r; i++)
  {
    printf("[ ");
    for(j = 0; j < c; j++)
    {
      printf("%5d ", *A++);
    }
    printf("]\n");
  }
  printf("\n");
}

void multiply()
// mat['!'] = mat['!'] times mat['@'] modulo MOD
{
  int *X;
  int *Xij;
  int i;
  int j;
  int x;
  int r;
  int mid;
  int c;

  if (mat['@'] == NULL)
  {
    die("tried to multiply by a null pointer.");
  }

  if (mat['!'] == NULL)
  {
    mat['!'] = mat['@'];
    bif_rows['!'] = bif_rows['@'];
    bif_cols['!'] = bif_cols['@'];
    mat['@'] = NULL;
    return;
  }

  if (bif_cols['!'] != bif_rows['@'])
  {
    die("tried to multiply matrices with mismatched dimensions.");
  }

  r = fib[bif_rows['!']];
  mid = fib[bif_cols['!']];
  c = fib[bif_cols['@']];
  X = malloc(r * c * sizeof(int));
  Xij = X;
  for(i = 0; i < r; i++)
  {
    for(j = 0; j < c; j++)
    {
      *Xij = 0;
      for(x = 0; x < mid; x++)
      {
        // add mat['!'][i,x]*mat['@'][x,j]
        *Xij += mat['!'][i*mid + x] * mat['@'][x*c + j];
	      *Xij %= MOD;
      }
      Xij++;
    }
  }
  free(mat['!']);
  mat['!'] = X;
  bif_cols['!'] = bif_cols['@'];
}

void tensor(int y[], int bif_y_r, int bif_y_c)
// mat['@'] = mat['@'] fibonacci-tensor y
{
  int *X;
  int *Xij;
  int iX;
  int iY;
  int jX;
  int jY;
  int bif_r;
  int bif_c;

  if(mat['@'] == NULL)
  {
    mat['@'] = copy_array(y, fib[bif_y_r] * fib[bif_y_c]);
    bif_rows['@'] = bif_y_r;
    bif_cols['@'] = bif_y_c;
    return;
  }

  bif_r = bif_rows['@'] + bif_y_r - 3;
  bif_c = bif_cols['@'] + bif_y_c - 3;
  X = malloc(fib[bif_r] * fib[bif_c] * sizeof(int));

  Xij = X;
  for(iX = 0; iX < fib[bif_rows['@']]; iX++)
  {
    for(iY = (fibword[iX] ? fib[bif_y_r - 1] : 0);
        iY < (fibword[iX] ? fib[bif_y_r] : fib[bif_y_r - 1]);
	      iY++)
    {
      for(jX = 0; jX < fib[bif_cols['@']]; jX++)
      {
        for(jY = (fibword[jX] ? fib[bif_y_c - 1] : 0);
            jY < (fibword[jX] ? fib[bif_y_c] : fib[bif_y_c - 1]);
            jY++)
        {
          *Xij = mat['@'][iX * fib[bif_cols['@']] + jX] * y[iY * fib[bif_y_c] + jY];
          *Xij %= MOD;
          Xij++;
        }
      }
    }
  }
  bif_rows['@'] = bif_r;
  bif_cols['@'] = bif_c;
  free(mat['@']);
  mat['@'] = X;
}

void exec_char(char c)
{
  int bif_r;
  int bif_c;

  if (c == ' ')  //ignore spaces
  {
    return;
  }

  if (c == '\n')
  {
    if ((mat['!'] == NULL) && (mat['@'] == NULL)) // 3rd carriage return
    {
      return;
    }
    if (mat['@'] == NULL) // 2nd carriage return
    {
      mat['!'] = NULL;
      return;
    }
    if (mat['!'] == NULL)
    {
      mat['!'] = mat['@'];
	    bif_rows['!'] = bif_rows['@'];
	    bif_cols['!'] = bif_cols['@'];
      mat['@'] = NULL;
      return;
	  }
	  multiply();
	  mat['@'] = NULL;
	  return;
  }

  if(('A' <= c) && (c <= 'Z'))
  {
    c += 'a' - 'A';
    bif_r = bif_rows['!'];
    bif_c = bif_cols['!'];
    mat[c] = copy_array(mat['!'], fib[bif_r] * fib[bif_c]);
    bif_rows[c] = bif_r;
    bif_cols[c] = bif_c;
    return;
  }

  if(mat[c] == NULL)
  {
    die("unknown character.");
  }

  tensor(mat[c], bif_rows[c], bif_cols[c]);
}

int main()
{
  char tangle[MAXSTRING];
  char line[MAXSTRING];

  int i;

  initialize();
  tangle[0]='\0';
  line[0]='\0';
  while(line[0] != '.')
  {
    fgets(line, MAXSTRING, stdin);
    strcat(tangle,line);
  }

  for(i = 0; tangle[i] != '.'; i++)
  {
    exec_char(tangle[i]);
  }
  prettyprint();
  freeall();
}
 
