#pragma once
#include <string>
#include <vector>
#include <io.h> 

using namespace std;
class GetFiles
{
public:
	const vector<string>& getFiles(const string& path);
	//char *getcwd( char *buffer, int maxlen ��
	//�����ܹ���ȡ��ǰ�Ĺ���Ŀ¼��������˵�����Ὣ��ǰ����Ŀ¼�ľ���·�����Ƶ�����buffer��ָ���ڴ�ռ���,����maxlenΪbuffer�Ŀռ��С��
	string GetWorkingDirectory();
private:
	vector<string> m_FilesVector;
};

