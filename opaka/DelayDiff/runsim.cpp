#include <admodel.h>
#include <stdio.h>

int main(int argc, char const *argv[])
{
	int i;
	int n =100;
	adstring arg;
	for (int i = 0; i < n; ++i)
	{
			arg = "./ddmod -sim " +str(i) +" -nox -est";
			system(arg);
	}
	return 0;
}