#pragma once
#include <string>
#include <vector>
#include <io.h> 

using namespace std;
class GetFiles
{
public:
	const vector<string>& getFiles(const string& path);
	//char *getcwd( char *buffer, int maxlen ）
	//函数能够获取当前的工作目录，具体来说，它会将当前工作目录的绝对路径复制到参数buffer所指的内存空间中,参数maxlen为buffer的空间大小。
	string GetWorkingDirectory();
private:
	vector<string> m_FilesVector;
};

