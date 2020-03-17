#ifndef __GILLESPIE_H__
#define __GILLESPIE_H__
#include "kmc.h"
using namespace std;
class Gillespie
{
private:
	void DeleteCurrentItemLib();
	void DeletecLibrary();
	void DeletehLibrary();
	void DeleteAllChain();
	void DeleteItemConv();
	void DeleteReactantName();
	bool Initialization();
	bool InitcLibrary();
	void InithLibrary();
	void InitItemConv();
	void InitReactantName();
	void InitAllChainArray();
	void UpdatehLibrary();
private:
	void FindNumbersAndStoichiometric(const Reaction* const pCurrentReaction, const string& name, int& numbers, int& stoichiometric);
	void ChainPropagation(const Reaction* const pCurrentReaction);
	void DoPropagation(const string& strAddUnit, const string& strEndUnit, const string& strInitiator);
	double Cal_c_GMC(const double& k, const double& x, const double& y);
	double CalCurrentH(const Reaction* const pCurrentReaction);
	double GenerateT();
	int GenerateU();
	void ExecuteAReact(const Reaction* const pCurrentReaction);
	double Cal_Sum_hMultipyC(const unsigned int pos);
	void TerminateChainPropagation();
	bool IsMeetTermination();
	double CalculateConv(const string& strName);
	void ExecuteStepAdder();
	void StepForward();
	void  TakeSample(double& time_pre);
	vector<string>* GetReactantName(const Reaction* const pCurrentReaction, vector<string>* pReactantNameArray = NULL);
	void WriteConvToFile();
public:
	void StartSimulation();
	vector<pPerChain>* GetAllChain() { 
		return m_pAllChain; 
	}
    Gillespie(KmcDoc* pKmcDoc = NULL);
    ~Gillespie();
private:
	vector<pPerChain>* m_pAllChain;
	vector<SubstanceItem*>* m_pCurrentItemLib;
	vector<double>* m_pcLibrary;
	vector<double>* m_phLibrary;
	deque<double>* m_pItemConv;
	vector<string>* m_pReactantName;
	KmcDoc* m_pKmcDoc;
	int m_StepsAddXuhao;
	double m_CurrentTime;
	vector<SpecieVector*> m_Specie;
};
#endif