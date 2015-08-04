GLOBALS_SECTION
 #include <andre.inc>

DATA_SECTION
 int N;
 
 vector yy1(1,5);
 !! yy1(1) = 1; yy1(2) = 1; yy1(3) = 1; yy1(4) = 1; yy1(5) = 10;
 
PARAMETER_SECTION 
 init_number XX;
 vector yy2(1,5);
 !! yy2 = yy1;
 
 objective_function_value obj_fun;
 
PROCEDURE_SECTION
 cout << sum(square(yy1)) << endl;
 cout << square2(yy1,5.0) << endl;
 cout << square2(yy2,5.0) << endl;
 cout << square2(yy2,5.0,5) << endl;

 obj_fun = XX*XX;
 
 
