#include "analysistool.h"
#include "kmcdoc.h"
#include "reaction.h"
#include <fstream>
using namespace std;

AnalysisTool::AnalysisTool(vector<pPerChain>* pAllChain, KmcDoc* pKmcDoc) :
	m_cal_Mn(0),
	m_cal_Mw(0),
	m_cal_PDI(0),
	m_pfAllBlock(NULL),
	m_pKmcDoc(pKmcDoc),
	m_pAllChain(pAllChain),
	m_pBlockNames(NULL),
	m_pMnList(NULL),
	m_pfBlocksInfo(NULL)
{ 
}

void AnalysisTool::SetAllChainandKmcDoc(vector<pPerChain>* pAllChain, KmcDoc* pKmcDoc)
{
	m_pKmcDoc = pKmcDoc;
	m_pAllChain = pAllChain;
}

void AnalysisTool::DeleteBlocksInfo()
{
	if (NULL == m_pfBlocksInfo)
		return;
	for (auto it = m_pfBlocksInfo->begin(); it != m_pfBlocksInfo->end(); it++)
	{
		for (auto it2 = it->_pfBlocks->begin(); it2 != it->_pfBlocks->end(); it2++)
		{
			for (auto it3 = (*it2)->begin(); it3 != (*it2)->end(); it3++)
			{
				delete (*it3);
			}
			delete (*it2);
		}
		delete it->_pfBlocks;
	}
	delete m_pfBlocksInfo;
	m_pfBlocksInfo = NULL;
}

void AnalysisTool::DeleteMnList()
{
	if (NULL == m_pMnList)
		return;
	delete m_pMnList;
	m_pMnList = NULL;
}

void AnalysisTool::DeleteBlockName()
{
	if (NULL == m_pBlockNames)
		return;
	delete m_pBlockNames;
	m_pBlockNames = NULL;
}

void AnalysisTool::Delete_f_AllBlock()
{
	if (NULL == m_pfAllBlock)
		return;
	for (auto itor = m_pfAllBlock->begin(); itor != m_pfAllBlock->end(); itor++)
	{
		size_t ncount = (*itor)->size();
		for (size_t i = 0; i < ncount; i++)
		{
			if(NULL != (*itor)->at(i)->pLengthArray)
				delete (*itor)->at(i)->pLengthArray;
			if (NULL != (*itor)->at(i)->pPositionArray)
				delete (*itor)->at(i)->pPositionArray;
			if (NULL != (*itor)->at(i)->pPositionDistributionArray)
				delete (*itor)->at(i)->pPositionDistributionArray;
			if (NULL != (*itor)->at(i)->pLengthDistributionArray)
				delete (*itor)->at(i)->pLengthDistributionArray;
			delete (*itor)->at(i);
		}
		delete *itor;
	}
	delete m_pfAllBlock;
	m_pfAllBlock = NULL;
}

AnalysisTool::~AnalysisTool()
{
	DeleteBlockName();
	Delete_f_AllBlock();
	DeleteMnList();
	DeleteBlocksInfo();
}
bool AnalysisTool::Init()
{
	if (NULL == m_pAllChain)
		return false;
	Delete_f_AllBlock();
	m_pfAllBlock = new vector<pBlockFloatArray>;
	DeleteMnList();
	m_pMnList = new vector<Distribution>;
	InitBlock_Name_Pair_Info();
	return true;
}
void AnalysisTool::InitBlock_Name_Pair_Info()
{
	DeleteBlockName();
	DeleteBlocksInfo();
	m_pBlockNames = new vector<string>;
	m_pfBlocksInfo = new vector<fBlockInfo>;
	for (auto itor = m_pKmcDoc->GetReactionSet().begin(); itor != m_pKmcDoc->GetReactionSet().end(); itor++)
	{		
		string	strAddname = (*itor)->GetReactionInfo()._strAddUnit;
		if (strAddname.empty())
			continue;
		bool nFlag = false;
		for (auto itBlock = m_pBlockNames->begin(); itBlock != m_pBlockNames->end(); itBlock++)
		{
			if (strAddname == *itBlock){
				nFlag = true;
				break;
			}
		}
		if (false == nFlag)
		{
			m_pBlockNames->push_back(strAddname);
			m_BlkNamePair.push_back(make_pair(strAddname, number_with_sameName(0)));
			m_pfBlocksInfo->push_back(fBlockInfo(strAddname, max_count(0), new vector<vector<BlockItem*>*>));
		}
	}
}

void AnalysisTool::AnalyzeAChain(pPerChain pChain)
{
	string strNameKeep = pChain->at(0);
	int iLength = 0;
	int iPosInChain = 1;
	for (auto it_Info = m_pfBlocksInfo->begin(); it_Info != m_pfBlocksInfo->end(); it_Info++) {
		(*it_Info)._pfBlocks->push_back(new vector<BlockItem*>);
	}
	for (auto it_pair = m_BlkNamePair.begin(); it_pair != m_BlkNamePair.end(); it_pair++) {
		(*it_pair).second = 0;	//初始化 m_BlkNamePair 中的 iNumberWithSameName 为 0.
	}	
	for (auto itor = pChain->begin(); itor != pChain->end(); itor++)
	{
		if (strNameKeep == *itor)
			iLength++;
		else{
			size_t iNumberWithSameName = 0;
			for (auto it_pair = m_BlkNamePair.begin(); it_pair != m_BlkNamePair.end(); it_pair++){
				if (strNameKeep == (*it_pair).first)
					iNumberWithSameName = (*it_pair).second++;
			}
			for (auto it_Info = m_pfBlocksInfo->begin(); it_Info != m_pfBlocksInfo->end(); it_Info++){
				if (strNameKeep == (*it_Info)._name)
				{
					if ((*it_Info)._maxcount < (iNumberWithSameName + size_t(1)))
						(*it_Info)._maxcount = iNumberWithSameName + size_t(1);
					(*it_Info)._pfBlocks->back()->push_back(new BlockItem(strNameKeep, iLength, iPosInChain, iNumberWithSameName));
				}
			}
			strNameKeep = *itor;
			iPosInChain = (int)(itor - pChain->begin()) + 1;
			iLength = 1;
		}
	}
}

pBlockFloatArray AnalysisTool::EvaluateSequence(const fBlockInfo& blockInfo)
{
	pBlockFloatArray pBlockList = new vector<BlockFloat*>;
	for (size_t iNum = 0; iNum < blockInfo._maxcount; iNum++)
	{
		BlockFloat* pfBlock = new BlockFloat;
		pfBlock->blockName = blockInfo._name;
		pfBlock->pLengthDistributionArray = new vector<Distribution>;
		pfBlock->pLengthDistributionArray->reserve(blockInfo._pfBlocks->size());
		pfBlock->pPositionDistributionArray = new vector<Distribution>;
		pfBlock->pPositionDistributionArray->reserve(blockInfo._pfBlocks->size());
		pfBlock->pLengthArray = new vector<double>;
		pfBlock->pLengthArray->reserve(blockInfo._pfBlocks->size());
		pfBlock->pPositionArray = new vector<double>;
		pfBlock->pPositionArray->reserve(blockInfo._pfBlocks->size());
		for (auto it_Chain = blockInfo._pfBlocks->begin(); it_Chain != blockInfo._pfBlocks->end(); it_Chain++){
			if (iNum > (*it_Chain)->size() || iNum == (*it_Chain)->size())
				continue;
			pfBlock->blockNumberWithSameName = (*it_Chain)->at(iNum)->blockNumberWithSameName;
			pfBlock->valueableChainNumber++;
			pfBlock->blockLength += (double)(*it_Chain)->at(iNum)->blockLength;
			pfBlock->pLengthArray->push_back((double)(*it_Chain)->at(iNum)->blockLength);
			pfBlock->blockPosInChain += (double)(*it_Chain)->at(iNum)->blockPosInChain;
			pfBlock->pPositionArray->push_back((double)(*it_Chain)->at(iNum)->blockPosInChain);
			JudgeDistribution((double)(*it_Chain)->at(iNum)->blockLength, pfBlock->pLengthDistributionArray);
			JudgeDistribution((double)(*it_Chain)->at(iNum)->blockPosInChain, pfBlock->pPositionDistributionArray);
		}
		pfBlock->blockLength /= (double)pfBlock->valueableChainNumber;
		pfBlock->blockPosInChain /= (double)pfBlock->valueableChainNumber;
		double sumLen = 0, sumPos = 0;
		size_t size = pfBlock->pLengthArray->size();
		for (size_t i = 0; i < size; i++){
			sumLen += pow(pfBlock->pLengthArray->at(i) - pfBlock->blockLength,2);
			sumPos += pow(pfBlock->pPositionArray->at(i) - pfBlock->blockPosInChain, 2);
		}
		pfBlock->lengthStandardDeviation = sqrt(sumLen / pfBlock->valueableChainNumber);
		pfBlock->positionStandardDeviation = sqrt(sumPos / pfBlock->valueableChainNumber);
		pBlockList->push_back(pfBlock);
	}
	return pBlockList;
}
void AnalysisTool::JudgeDistribution(const double value, vector<Distribution>* pDistributionArray)
{
	int nFlag = 0;
	for(auto itor = pDistributionArray->begin(); itor != pDistributionArray->end(); itor++){
		nFlag = 0;
		if(value == (*itor).fValue){
			nFlag = 1;
			(*itor).iNumbers++;
			break;
		}
	}
	if (!nFlag)
		pDistributionArray->push_back(Distribution(value,1));
}

double AnalysisTool::GetSubtanceItemWeight(const string& name)
{
	for (auto itor = m_pKmcDoc->GetSubtanceItemLibrary()->begin(); itor != m_pKmcDoc->GetSubtanceItemLibrary()->end(); itor++)
	{
		if (name == (*itor)->substanceName)
			return (*itor)->substanceMolecularWeight;
	}
	return 0;
}

void AnalysisTool::AnalyzeMnandPDI()
{
	vector<double>* pMn = new vector<double>;
	for (auto itor = m_pAllChain->begin(); itor != m_pAllChain->end(); itor++)	{
		double fMn = 0;
		for (auto itor2 = (*itor)->begin(); itor2 != (*itor)->end(); itor2++)
			fMn += GetSubtanceItemWeight((*itor2));
		fMn = (double)((size_t)(fMn / 100) * 100);
		pMn->push_back(fMn);
	}
	for (auto itor = pMn->begin(); itor != pMn->end(); itor++)
		JudgeDistribution((*itor), m_pMnList);
	delete pMn;
	double Mn_Num(0), Mn_Mn_Num(0), num_Mn_Num(0);
	m_cal_PDI = m_cal_Mw = m_cal_Mn = 0;
	for (auto itor = m_pMnList->begin(); itor != m_pMnList->end(); itor++){
		Mn_Num += (*itor).fValue * (*itor).iNumbers;
		num_Mn_Num += (double)(*itor).iNumbers;
	}
	for (auto itor = m_pMnList->begin(); itor != m_pMnList->end(); itor++){
		Mn_Mn_Num = (*itor).fValue * (*itor).fValue * (double)(*itor).iNumbers;
		if (Mn_Num)
			m_cal_Mw += (Mn_Mn_Num / (double)Mn_Num);
	}
	if (num_Mn_Num > 0)
		m_cal_Mn = Mn_Num / num_Mn_Num;
	if (m_cal_Mn)
		m_cal_PDI = m_cal_Mw / m_cal_Mn;
	cout << "\t[ Analysis Mn and PDI ] done" << endl;
}

void AnalysisTool::AnalyzeSequence()
{
	for (auto itor = m_pAllChain->begin(); itor != m_pAllChain->end(); itor++)
		if (!(*itor)->empty())
			AnalyzeAChain((*itor));
	for (auto itor2 = m_pfBlocksInfo->begin(); itor2 != m_pfBlocksInfo->end(); itor2++)
		m_pfAllBlock->push_back(EvaluateSequence(*itor2));
	cout << "\t[ Analysis Sequence ] done" << endl;
}

void AnalysisTool::StartAnalysis()
{
	if (false == Init())
		return;
	cout << endl << "[ Analysis ] is running..." << endl << endl;
	AnalyzeSequence();
	AnalyzeMnandPDI();
	cout << endl << "[ Analysis ] done!" << endl << endl;
	WriteResultToFile();
}

void AnalysisTool::WriteResultToFile()
{
	ofstream outFile;//创建了一个ofstream 对象
	string stringOutFile = m_pKmcDoc->GetSimuFileName() + ".out";
	outFile.open(stringOutFile.data(), ios::app);//outFile 与一个文本文件关联
	outFile << endl;
	for (auto itor = m_pfAllBlock->begin(); itor != m_pfAllBlock->end(); itor++)
	{
		for (auto itor2 = (*itor)->begin(); itor2 != (*itor)->end(); itor2++)
		{
			outFile << "Block: " << (*itor2)->blockName << " \t#" << (*itor2)->blockNumberWithSameName + 1 << ":" << endl;
			outFile << "Average Length: " << (*itor2)->blockLength << endl << "Average Position: " << (*itor2)->blockPosInChain << endl;
			outFile << "Length Standard Deviation: " << (*itor2)->lengthStandardDeviation << endl << "Position Standard Deviation: " << (*itor2)->positionStandardDeviation << endl;
			outFile << "Length Distribution:" << endl << "Length" << "\tCount" << endl;
			for (auto itor3 = (*itor2)->pLengthDistributionArray->begin(); itor3 != (*itor2)->pLengthDistributionArray->end(); itor3++)
				outFile << (*itor3).fValue << "\t" << (*itor3).iNumbers << endl;
			outFile << endl;
			outFile << "Position Distribution:" << endl << "Position" << "\tCount" << endl;
			for (auto itor3 = (*itor2)->pPositionDistributionArray->begin(); itor3 != (*itor2)->pPositionDistributionArray->end(); itor3++)
				outFile << (*itor3).fValue << "\t" << (*itor3).iNumbers << endl;
			outFile << endl;
		}
	}
	outFile << endl << "Molecular weight distribution: " << endl << "Weight" << "\tCount" << endl;
	for (auto itor = m_pMnList->begin(); itor != m_pMnList->end(); itor++)
		outFile << (*itor).fValue << "\t" << (*itor).iNumbers << endl;
	outFile << endl << "Mn = " << m_cal_Mn << endl << "Mw = " << m_cal_Mw << endl << "PDI = " << m_cal_PDI << endl;
	outFile.close();
}
