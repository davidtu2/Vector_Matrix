Name: David Tu
Contact: david.tu2@csu.fullerton.edu / 626-497-3531

Please refer to the pdf file for source code listing as well as output results

Program Description: This is a program which tests 3D Matricies and Vectors by using assertions. The program will check for correct mathematical computations for
the matricies and vectors. If they are all correct, then the program will indicate that all assertions have passed and exists without error

How to Run the Program:
	Pre-Requisite(s): Installation of Visual Studio 2015
	1. Open Visual Studio and create an empty C++ project (You may need to install C++ if this is your first time)
	2. Copy/Paste the matrix3d_t.h, vector3d_t.h and main.cpp files into the newly crated empty project directory
	3. On the Solution Explorer window in Visual Studio, right click Header Files > Add > Existing Item to add
	the vector3d_t.h and matrix3d_t.h files into your solution
	4. On the Solution Explorer window in Visual Studio, right click Cource Files > Add > Existing Item to add
	main.cpp into your solution
	5. Go to Build > Build Solution
	6. Then click on Local Windows Debugger (the "Play" button) to run the program

Known Bugs/Issues: If an assertion fails, the program will error and terminate due to no exception handling mechanism as of now