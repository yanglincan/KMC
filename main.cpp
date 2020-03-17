//#include "kmc.h"
//#include "reaction.h"
//#include "kmcfile.h"
#include "kmcdoc.h"
#include "GetFiles.h"
//#include <thread>
//#include <time.h>
//#include <fstream>

//#ifdef _DEBUG  
//#define New   new(_NORMAL_BLOCK, __FILE__, __LINE__)  
//#endif  
//
//#define CRTDBG_MAP_ALLOC    
//#include <stdlib.h>    
//#include <crtdbg.h>  

using namespace std;

const string GetSimuFile()
{
	cout << endl << "-> Place the simulation file into the foldor that stores kmc program. <-" << endl << endl;
	GetFiles filename;
	const vector<string>& v = filename.getFiles(filename.GetWorkingDirectory());
	for_each(v.begin(), v.end(), [](string filename){static int i = 1; cout << "#" << i++ << ": " << filename << endl;	});
	cout << endl << "-> Input and select [e.g. ""1 + enter""] <-" << endl << endl;
	size_t select = 1;
	cin >> select;
	return v.at(select - 1);
}

int main()
{
	{
		KmcDoc doc(GetSimuFile());

		/*std::random_device rd;
		cout << rd.min() << "  " << rd.max() << endl;
		for (int i = 0; i < 100; i++)
			cout << rd()%10+1 << " ";
		cout << endl;*/
	}
	//_CrtDumpMemoryLeaks();
	system("pause");
	return 0;
}
