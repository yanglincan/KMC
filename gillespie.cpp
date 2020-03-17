#include "gillespie.h"
#include "kmcdoc.h"
#include "reaction.h"
#include <time.h>
#include "analysistool.h"
#include <fstream>
#include <random>
using namespace std;

int CreateRandomInteger(int nMax)
{
	if (nMax < 0 || nMax > 300000)
		return -1;
	if (0 == nMax)
		return 0;
	return ((rand() * 10) + (rand() % 10)) % nMax + 1;	//取得[a,b]的随机整数：rand()%(b-a+1)+a, RAND_MAX = 32767
}

double CreateRandomDouble()
{
	return (rand() / (RAND_MAX + 1.0));
}

Gillespie::Gillespie(KmcDoc* pKmcDoc) :
	m_pAllChain(NULL),
	m_pCurrentItemLib(NULL),
	m_pcLibrary(NULL),
	m_phLibrary(NULL),
	m_pItemConv(NULL),
	m_pReactantName(NULL),
	m_pKmcDoc(pKmcDoc),
	m_StepsAddXuhao(1),
	m_CurrentTime(0)
{
}

Gillespie::~Gillespie()
{
	DeletecLibrary();
	DeletehLibrary();
	DeleteCurrentItemLib();
	DeleteAllChain();
	DeleteItemConv();
	DeleteReactantName();
	for (auto it = m_Specie.begin(); it != m_Specie.end(); it++)
		delete (*it);
}

void Gillespie::DeleteCurrentItemLib()
{
	if (NULL == m_pCurrentItemLib)
		return;
	for (auto itor = m_pCurrentItemLib->begin(); itor != m_pCurrentItemLib->end(); itor++)
		delete *itor;
	delete m_pCurrentItemLib;
	m_pCurrentItemLib = NULL;
}

void Gillespie::DeletecLibrary()
{
	if (NULL == m_pcLibrary)
		return;
	delete m_pcLibrary;
	m_pcLibrary = NULL;
}

void Gillespie::DeletehLibrary()
{
	if (NULL == m_phLibrary)
		return;
	delete m_phLibrary;
	m_phLibrary = NULL;
}

void Gillespie::DeleteAllChain()
{
	if (NULL == m_pAllChain)
		return;
	for (auto itor = m_pAllChain->begin(); itor != m_pAllChain->end(); itor++)
		delete *itor;
	delete 	m_pAllChain;
	m_pAllChain = NULL;
}

void Gillespie::DeleteReactantName()
{
	if (NULL == m_pReactantName)
		return;
	delete m_pReactantName;
	m_pReactantName = NULL;
}

void Gillespie::DeleteItemConv()
{
	if (NULL == m_pItemConv)
		return;
	delete 	m_pItemConv;
	m_pItemConv = NULL;
}

bool Gillespie::Initialization()
{
	m_StepsAddXuhao = 1;
	InitReactantName();
	InithLibrary();
	InitItemConv();
	if (false == InitcLibrary())
		return false;
	InitAllChainArray();
	return true;
}

bool Gillespie::InitcLibrary()
{
	DeletecLibrary();
	m_pcLibrary = new  vector<double>;
	for (auto itor = m_pKmcDoc->GetReactionSet().begin(); itor != m_pKmcDoc->GetReactionSet().end(); itor++)
	{//_finite(t)返回1，未溢出；返回0，溢出
		double cValue = Cal_c_GMC(
			(*itor)->GetReactionRateConstant(),
			(*itor)->GetRelationBetween_kMC_k().x,
			(*itor)->GetRelationBetween_kMC_k().y);
		if (_finite(cValue))
			m_pcLibrary->push_back(cValue);
		else
			return false;
	}
	return true;
}

void Gillespie::InitReactantName()
{
	DeleteReactantName();
	m_pReactantName = new vector<string>;
	for (auto itor = m_pKmcDoc->GetReactionSet().begin(); itor != m_pKmcDoc->GetReactionSet().end(); itor++)
		GetReactantName((*itor), m_pReactantName);
}

void Gillespie::InithLibrary()
{
	DeleteCurrentItemLib();
	m_pCurrentItemLib = new vector<SubstanceItem*>;
	for (auto itor = m_pKmcDoc->GetSubtanceItemLibrary()->begin(); itor != m_pKmcDoc->GetSubtanceItemLibrary()->end(); itor++)
		m_pCurrentItemLib->push_back(new SubstanceItem(
			(*itor)->substanceName,
			(*itor)->substanceAmount,
			(*itor)->substanceisReactant,
			(*itor)->substanceisProduct,
			(*itor)->substanceMolecularWeight
		));
	DeletehLibrary();
	m_phLibrary = new  vector<double>;
	for (auto itor = m_pKmcDoc->GetReactionSet().begin(); itor != m_pKmcDoc->GetReactionSet().end(); itor++)
		m_phLibrary->push_back(CalCurrentH(*itor));
}

void Gillespie::InitItemConv()
{
	DeleteItemConv();
	m_pItemConv = new  deque<double>;
}

void Gillespie::UpdatehLibrary(void)
{
	size_t size = m_phLibrary->size();
	for (size_t i = 0; i < size; i++)
		m_phLibrary->at(i) = CalCurrentH(m_pKmcDoc->GetReactionSet().at(i));
}

void Gillespie::InitAllChainArray()
{
	DeleteAllChain();
	string strInitName, strSpecie;
	for (auto itor = m_pKmcDoc->GetReactionSet().begin(); itor != m_pKmcDoc->GetReactionSet().end(); itor++)
	{
		if (!(*itor)->GetReactionInfo()._strInitiator.empty())
			strInitName = strSpecie =(*itor)->GetReactionInfo()._strInitiator;
		if (!(*itor)->GetReactionInfo()._strInitiatSpecie.empty())
			strSpecie = (*itor)->GetReactionInfo()._strInitiatSpecie;
		if (!strSpecie.empty())
		{
			bool nFlag = false;
			for(auto itSpecie = m_Specie.begin(); itSpecie != m_Specie.end(); itSpecie++)
				if(strSpecie == (*itSpecie)->_SpecieName)
					nFlag = true;
			if (false == nFlag)
				m_Specie.push_back(new SpecieVector(strSpecie, new vector<int>));
		}
	}
	int nSize = 0;
	for (auto itor3 = m_pKmcDoc->GetSubtanceItemLibrary()->begin(); itor3 != m_pKmcDoc->GetSubtanceItemLibrary()->end(); itor3++)
	{
		if (strInitName == (*itor3)->substanceName)
		{
			nSize = (*itor3)->substanceAmount;
			m_pAllChain = new vector<pPerChain>;
			for (auto itSpecie = m_Specie.begin(); itSpecie != m_Specie.end(); itSpecie++)
				(*itSpecie)->_pVector->reserve(nSize);
			break;
		}
	}
	for (int i = 0; i < nSize; i++)
		m_pAllChain->push_back(new vector<string>);

	for (auto itSpecie = m_Specie.begin(); itSpecie != m_Specie.end(); itSpecie++)
	{
		if (strInitName == (*itSpecie)->_SpecieName)
		{
			for (int i = 0; i < nSize; i++)
				(*itSpecie)->_pVector->push_back(i);
		}
	}
}

vector<string>* Gillespie::GetReactantName(const Reaction* const pCurrentReaction, vector<string>* pReactantNameArray)
{
	if(NULL == pReactantNameArray)
		pReactantNameArray = new vector<string>;
	//调用该函数，要手动delete掉pReactantNameArray
	int nExistFlag = 0;
	for (auto itor = pCurrentReaction->GetElementVector().begin(); itor != pCurrentReaction->GetElementVector().end(); itor++)
	{
		if (REACTANT == (*itor)->elementType)
		{
			for (auto itor2 = pReactantNameArray->begin(); itor2 != pReactantNameArray->end(); itor2++)
			{
				if ((*itor)->elementName == *itor2)
				{
					nExistFlag = 1; 
					break;
				}
			}
			if (!nExistFlag)
			{
				pReactantNameArray->push_back((*itor)->elementName);
				nExistFlag = 0;
			}
		}
	}
	return pReactantNameArray;
}

double Gillespie::CalCurrentH(const Reaction* const pCurrentReaction)
{
	double finalValue = 1;	//finalValue即是h值
	vector<string>* pReactantNameArray = NULL;
	pReactantNameArray = GetReactantName(pCurrentReaction);
	for (auto itor = pReactantNameArray->begin(); itor != pReactantNameArray->end(); itor++)
	{
		int numbers = 0;
		int stoichiometric = 0;
		FindNumbersAndStoichiometric(pCurrentReaction, *itor, numbers, stoichiometric);
		if (1 == stoichiometric)
			numbers = numbers;
		else if (2 == stoichiometric)
			numbers = numbers * (numbers - 1) / 2;
		else if (3 == stoichiometric)
			numbers = numbers * (numbers - 1) * (numbers - 2) / 6;
		finalValue = finalValue * numbers;
	}
	if (NULL != pReactantNameArray)
		delete	pReactantNameArray;
	return finalValue > 0 ? finalValue : 0;
}

void Gillespie::FindNumbersAndStoichiometric(const Reaction* const pCurrentReaction, const string& name, int& numbers, int& stoichiometric)
{
	stoichiometric = 0;
	numbers = 0;
	for (auto itor = m_pCurrentItemLib->begin(); itor != m_pCurrentItemLib->end(); itor++)
	{
		if (name == (*itor)->substanceName)
		{ 
			numbers = (*itor)->substanceAmount; 
			break; 
		}
	}
	for (auto itor = pCurrentReaction->GetElementVector().begin(); itor != pCurrentReaction->GetElementVector().end(); itor++)
	{
		if(name == (*itor)->elementName && REACTANT == (*itor)->elementType)
			stoichiometric = stoichiometric + (*itor)->stoichiometricNumber;
	}
}

double Gillespie::Cal_c_GMC(const double& k, const double& x, const double& y)
{
	double na = 6.02 * pow((double)10, (int)23);
	return x * k / (pow(na * m_pKmcDoc->GetVolume(), (double)y));	//(x,y): k_mc = x*k/[(NaV)^y]。
}

void Gillespie::StartSimulation()
{
	if (NULL == m_pKmcDoc || !Initialization())
		return;
	cout << endl << "[ Simulation ] is running..." << endl << endl;
	double time_pre = 0;
	srand((unsigned)time(0));
	while (false == IsMeetTermination())
	{
		ExecuteStepAdder();
		StepForward();
		TakeSample(time_pre);
	}
	TerminateChainPropagation();
	cout << endl<<"[ Simulation ] done!" << endl << endl;
	WriteConvToFile();
}

bool Gillespie::IsMeetTermination()
{
	return (CalculateConv(m_pKmcDoc->GetTerminatrName()) > m_pKmcDoc->GetTerminatrCove()) ? true : false;
}

void  Gillespie::TakeSample(double& time_pre)
{
	static long nSampleCount = 0;
	if ((m_CurrentTime - time_pre) > m_pKmcDoc->GetIntervalTime())
	{
		m_pItemConv->push_back(m_CurrentTime);
		for (auto itor = m_pReactantName->begin(); itor != m_pReactantName->end(); itor++)
			m_pItemConv->push_back(CalculateConv((*itor)));
		m_pItemConv->push_back(-2);	// 每一轮结尾放 -2 ，用于后续分行输出
		cout << "\tSampling: "<< ++nSampleCount << endl;
		time_pre = m_CurrentTime;
	}
}

double Gillespie::CalculateConv(const string& strName)
{
	size_t nCount = m_pCurrentItemLib->size();
	for (size_t i = 0; i < nCount; i++)
	{
		if (strName == m_pCurrentItemLib->at(i)->substanceName)
		{
			int n = m_pCurrentItemLib->at(i)->substanceAmount;
			int n_0 = m_pKmcDoc->GetSubtanceItemLibrary()->at(i)->substanceAmount;
			return (0 == n_0) ? -1 : (((double)n_0 - (double)n)) / (double)n_0;
		}
	}
	return -1;
}

void Gillespie::ExecuteStepAdder()
{
	int nFlag = 0;
	double currentconv = CalculateConv(m_pKmcDoc->GetTerminatrName());
	for (auto itor = m_pKmcDoc->GetStepsAdder().begin(); itor != m_pKmcDoc->GetStepsAdder().end(); itor++)
	{
		//从m_StepsAddXuhao=1的序号开始添加，然后判断是否已经添加，若没有，查看转化率是否合适，合适则添加，否则扫描下一条目
		if (m_StepsAddXuhao == (*itor)->stepsNumber && false == (*itor)->stepsDone && (*itor)->stepsConv < currentconv)
		{
			size_t nCount = m_pCurrentItemLib->size();
			for (size_t i = 0; i < nCount; i++)
			{
				if ((*itor)->stepsName == m_pCurrentItemLib->at(i)->substanceName)
				{
					m_pCurrentItemLib->at(i)->substanceAmount +=  (*itor)->stepsAmount;
					m_pKmcDoc->GetSubtanceItemLibrary()->at(i)->substanceAmount += (*itor)->stepsAmount;
					(*itor)->stepsDone = true;
					nFlag++;
				}
			}
		}
	}
	if (nFlag) 
		m_StepsAddXuhao++;	//如果同一序号下已经有一次添加，则进行下一序号的添加。因此，同一序号下，转化率较大的加料条目会被忽略。
}

double Gillespie::GenerateT()
{
	double a = Cal_Sum_hMultipyC((int)m_pKmcDoc->GetReactionSet().size());
	if (-1 == a)
		return -1;
	double t = 0;
	int reCount = 0;
	do
	{
		if (20 < reCount++)
			m_pKmcDoc->Exit("Error");
		t = (1 / a) * log((double)(1 / CreateRandomDouble()));
	} while (!_finite(t));
	return t;
}

int Gillespie::GenerateU()
{
	double random = 0;
	do
		random = CreateRandomDouble();
	while (0 == random || 1 == random);
	int pos = (int)m_pKmcDoc->GetReactionSet().size();
	double a = Cal_Sum_hMultipyC((int)pos);	//计算所有反应的h乘以c之和
	for (int i = 0; i < pos; i++)
	{
		double a_left = Cal_Sum_hMultipyC(i);
		double a_right = Cal_Sum_hMultipyC(i + 1);
		double ra = random * a;
		if (a_left < ra && ra < a_right || ra == a_right)
			return i;
	}
	return -1;
}

void Gillespie::ExecuteAReact(const Reaction* const pCurrentReaction)
{
	for (auto itor = pCurrentReaction->GetElementVector().begin(); itor != pCurrentReaction->GetElementVector().end(); itor++)
	{
		for (auto itor2 = m_pCurrentItemLib->begin(); itor2 != m_pCurrentItemLib->end(); itor2++)
		{
			if ((*itor)->elementName == (*itor2)->substanceName)
				(*itor2)->substanceAmount += (((REACTANT == (*itor)->elementType) ? -1 : 1) * (*itor)->stoichiometricNumber);
		}
	}
}

void Gillespie::StepForward()
{
	double t = GenerateT();
	int u = GenerateU();
	if (-1 == t || -1 == u)
		return;
	m_CurrentTime = m_CurrentTime + t;
	ExecuteAReact(m_pKmcDoc->GetReactionSet().at(u));
	ChainPropagation(m_pKmcDoc->GetReactionSet().at(u));
	UpdatehLibrary();
}

void Gillespie::DoPropagation(const string& strAddUnit, const  string& strEndUnit, const string& strSpecie)
{
	vector<string>* pChain = NULL;
	int random = -1;
	vector<SpecieVector*>::iterator selectedIter;
	for (auto itSpecie = m_Specie.begin(); itSpecie != m_Specie.end(); itSpecie++)
	{
		if (strSpecie == (*itSpecie)->_SpecieName)
		{
			selectedIter = itSpecie;
			random = CreateRandomInteger(int((*selectedIter)->_pVector->size() - 1));
			break;
		}
	}
	if (random < 0)
		return;
	pChain = m_pAllChain->at((*selectedIter)->_pVector->at(random));
	pChain->push_back(strAddUnit);
	::swap((*selectedIter)->_pVector->at(random), (*selectedIter)->_pVector->back());
	for (auto itSpecie = m_Specie.begin(); itSpecie != m_Specie.end(); itSpecie++)
	{
		if (strEndUnit == (*itSpecie)->_SpecieName)
		{
			(*itSpecie)->_pVector->push_back((*selectedIter)->_pVector->back());
			break;
		}
	}
	(*selectedIter)->_pVector->pop_back();
}

void Gillespie::ChainPropagation(const Reaction* const  pCurrentReaction)
{
	DoPropagation(
		pCurrentReaction->GetReactionInfo()._strAddUnit,
		pCurrentReaction->GetReactionInfo()._strEndUnit, 
		pCurrentReaction->GetReactionInfo()._strInitiator.empty() ?
		pCurrentReaction->GetReactionInfo()._strInitiatSpecie : pCurrentReaction->GetReactionInfo()._strInitiator);
}

void Gillespie::TerminateChainPropagation()
{
	if (!m_pAllChain)
		return;
	size_t nCountAllChain = m_pAllChain->size();
	for (size_t i = 0; i < nCountAllChain; i++)
	{
		vector<string>* pAChain = m_pAllChain->at(i);
		if (!pAChain->size())
			continue;
		pAChain->push_back("END");
	}
}

double Gillespie::Cal_Sum_hMultipyC(const unsigned int pos)
{
	if (m_pcLibrary->size() != m_phLibrary->size() || pos > m_phLibrary->size())
		return -1;
	double sum = 0;
	for (unsigned int i = 0; i < pos; i++)
		sum += m_pcLibrary->at(i) * m_phLibrary->at(i);
	return sum;
}

void Gillespie::WriteConvToFile()
{
	ofstream outFile;//创建了一个ofstream 对象
	string stringOutFile = m_pKmcDoc->GetSimuFileName() + ".out";
	outFile.open(stringOutFile.data(), ios::app);//outFile 与一个文本文件关联
	outFile << "time" << "\t";
	for (auto itor = m_pReactantName->begin(); itor != m_pReactantName->end(); itor++)
		outFile << (*itor) << "\t";
	outFile << endl;
	char tempChar[256] = { 0 };
	for (auto itor2 = m_pItemConv->begin(); itor2 != m_pItemConv->end(); itor2++)
		(-2 == *itor2)? outFile << endl : outFile << *itor2 << "\t";
	outFile.close();
}
