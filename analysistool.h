#ifndef __ANALYSISTOOL_H__
#define __ANALYSISTOOL_H__
#include "kmc.h"
using namespace std;

class AnalysisTool
{
public:
	AnalysisTool(vector<pPerChain>* pAllChain = NULL, KmcDoc* pKmcDoc = NULL);
	~AnalysisTool();
	void SetAllChainandKmcDoc(vector<pPerChain>* pAllChain, KmcDoc* pKmcDoc);
private:
	void DeleteBlockName();
	void Delete_f_AllBlock();
	void DeleteMnList();
	void DeleteBlocksInfo();
	bool Init();
	void InitBlock_Name_Pair_Info();
	void AnalyzeAChain(pPerChain pChain);
	pBlockFloatArray EvaluateSequence(const fBlockInfo& blockInfo);
	void JudgeDistribution(const double value, vector<Distribution>* pDistributionArray);
	void WriteResultToFile();
	double GetSubtanceItemWeight(const string& name);
	void AnalyzeSequence();
	void AnalyzeMnandPDI();
private:
	vector<pPerChain>* m_pAllChain;
	KmcDoc* m_pKmcDoc;
private:
	vector<pBlockFloatArray>* m_pfAllBlock;
	vector<string>* m_pBlockNames;
	vector<Distribution>* m_pMnList;
	vector<pair<string, number_with_sameName>> m_BlkNamePair;
	vector<fBlockInfo>* m_pfBlocksInfo;
	double m_cal_Mn;
	double m_cal_Mw;
	double m_cal_PDI;
public:
	void StartAnalysis();
};

#endif

