DATA_SECTION
  dll_int n
  dll_init_vector x(1,n)
  dll_init_vector y(1,n)

PARAMETER_SECTION
  dll_init_number a
  dll_init_number b
  dll_number output
  objective_function_value f

PROCEDURE_SECTION
  dvariable Pred;

  f = 0;
  for (int i=1; i<n; i++)
  {
   Pred = a + b*x(i);
   f += square(Pred-y(i));
  }
  output = f;
