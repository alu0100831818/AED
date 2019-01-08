

int& min(int &a, int &b)
{ return (a < b ? a : b);
}

int main()
{ int i = 1, j = 2;

  min(i, j) = 3; // ¿es esto posible?
  return 0;
}