#ifndef __KMC_H__
#define __KMC_H__
#include <string>
#include <vector>
#include <deque>
//#include <list>
#include <stack>
#include <iostream>
#include <stdexcept>
#include <ctime>
#include <algorithm>
//#include <functional>

using namespace std;
enum SubtanceType { REACTANT, PRODUCT };
enum SubtanceRole { UNSPECIFIED, MONOMER, SPECIE, INITIATOR };
enum SyntaxError { SYNTAX_PASS, NON_SMALL_BRACKET , NON_BRACKET, NON_BIG_BRACKET, NON_ARROW, NON_RATE, BRACKET_ERROR,
	NON_NAME, NON_STOICHIOMETRIC_NUMBER, NO_REACTION, MORE_THAN_THREE, ONLY_REACTANT_OR_PRODUCT};
enum ReactionOrder { _UNSPECIFIED_TYPE, _ZERO_TYPE, _A_TYPE, _AB_TYPE, _2A_TYPE, _ABC_TYPE, _A2B_YTPE, _3A_TYPE };

class Reaction;
class KmcDoc;
class Gillespie;
class AnalysisTool;
class KmcException;

struct ReactionInfo {
	string _strAddUnit; 
	string _strEndUnit;
	string _strInitiator;
	string _strInitiatSpecie;
	ReactionInfo(const string add, const  string end, const string initiator, const string specie):
		_strAddUnit(add),
		_strEndUnit(end),
		_strInitiator(initiator),
		_strInitiatSpecie(specie){}
	ReactionInfo() {}
};

struct SpecieVector
{
	string _SpecieName;
	vector<int>* _pVector;
	SpecieVector(string name = "", vector<int>* p = NULL) :_SpecieName(name), _pVector(p) {}
	~SpecieVector() {
		if (NULL != _pVector)
			delete _pVector;
	}
};

class KmcSyntaxException :public exception
{
public:
	char const* what() const
	{
		switch (m_ErrorCode)
		{
		case NO_REACTION: return "\n����Ӧ����ʽ���� ��Ӧ����ʽ  ��\n\n ʾ����\n\n(1)[A] + (1)[B] > (2)[C] {0.128}\n";
		case NON_SMALL_BRACKET: return "\n����Ӧ����ʽ���� ��ѧ����ϵ����ʶ ��( )�� ��\n";
		case NON_BRACKET: return "\n����Ӧ����ʽ���� ��Ӧ��or�������Ʊ�ʶ ��[ ]�� ��\n";
		case NON_BIG_BRACKET: return "\n����Ӧ����ʽ���� ��Ӧ���ʳ�����ʶ�� ��{ }�� ��\n";
		case NON_ARROW:  return "\n����Ӧ����ʽ���� ��Ӧָ���ʶ ��>�� ��\n";
		case NON_RATE: return "\n����Ӧ����ʽ����Ӧ���ʳ��� ��ֵ ����\n";
		case BRACKET_ERROR: return "\n����Ӧ����ʽ����( ) [ ] { }�����˳�� ����\n";
		case NON_NAME: return "\n����Ӧ����ʽ����Ӧ��or���� ���� Ϊ�գ�\n";
		case NON_STOICHIOMETRIC_NUMBER: return "\n����Ӧ����ʽ����ѧ����ϵ�� ��ֵ Ϊ�գ�\n";
		case MORE_THAN_THREE: return "\n����Ӧ����ʽ����֧�� ��Ӧ�ﻯѧ����ϵ��֮�� > 3 �ķ�Ӧ��\n";
		case ONLY_REACTANT_OR_PRODUCT: return "\n����Ӧ����ʽ��������� ��Ӧ�� or ���� ���䣡\n";
		default:
			break;
		}
		return NULL;
	}
	KmcSyntaxException(SyntaxError error = SYNTAX_PASS) :m_ErrorCode(error) {}
private:
	SyntaxError m_ErrorCode;
};
typedef struct ReactionElement
{
	string elementName;		//��������
	SubtanceType elementType;	//REACTANT or PRODUCT
	SubtanceRole elementRole;	//UNSPECIFIED, MONOMER, SPECIE or INITIATOR
	int stoichiometricNumber;	//��ѧ������
	ReactionElement(string name = "", SubtanceType type = REACTANT, SubtanceRole role = UNSPECIFIED, int number = 0):
		elementName(name),elementType(type),elementRole(role),stoichiometricNumber(number){}
} ReactionElement;

typedef struct OrderController
{
	//(x,y): k_mc = x*k/[(NaV)^y]��������Ӧ����
	double x;
	double y;
	OrderController(double xx = 0, double yy = 0) :x(xx), y(yy) {}
} OrderController;

typedef struct StepsAdder
{
	int stepsNumber;	//�ֲ����ϵ���ţ���1��ʼ��������
	double stepsConv;	//�ڵ�ǰת����ʱ��ʵʩ�ֲ�����
	string stepsName;	//�ֲ����ϵ���������
	int stepsAmount;	//�ֲ����ϵ�����
	bool stepsDone;		//�ֲ������Ƿ���ɵı��
	StepsAdder(string name, int number = 0, double conv = 0, int acount = 0, bool done = false) :
		stepsName(name), 
		stepsNumber(number),
		stepsConv(conv),
		stepsAmount(acount),
		stepsDone(done){}
	friend ostream& operator << (ostream& os, const StepsAdder& adder)
	{
		return os << "#" << adder.stepsNumber << ": Conv(" << adder.stepsConv << ")  "
			<<adder.stepsName << "  " << adder.stepsAmount << endl;
	}
} StepsAdder;

typedef struct SubstanceItem
{
	string substanceName;
	int substanceAmount;
	bool substanceisReactant;
	bool substanceisProduct;
	double substanceMolecularWeight;
	SubstanceItem(string name, int amount = 0, bool isReactant = false, bool isProduct = false, double weight = 0) :
		substanceName(name),
		substanceAmount(amount),
		substanceisReactant(isReactant),
		substanceisProduct(isProduct),
		substanceMolecularWeight(weight){}
	friend ostream& operator << (ostream& os, const SubstanceItem& item)
	{
		os << item.substanceName<<"  ";
		if (item.substanceisReactant) os << "  Reaction";
		if (item.substanceisProduct) os << "  Product";
		os << "  Weight:" << item.substanceMolecularWeight << "  Amount:" << item.substanceAmount << endl;
		return os;
	}
} SubstanceItem;

typedef struct BlockItem
{
	string blockName;
	int blockLength;
	int blockPosInChain;
	size_t blockNumberWithSameName;
	BlockItem(string name, int length, int posinchain, size_t numberwithsamename) :
		blockName(name), blockLength(length),blockPosInChain(posinchain),blockNumberWithSameName(numberwithsamename){}
} BlockItem;

typedef struct Distribution
{
	double fValue;
	int iNumbers;
	Distribution(double value, int num) : fValue(value), iNumbers(num) {};
} Distribution;

typedef struct BlockFloat
{
	string blockName;
	double blockLength;
	double blockPosInChain;
	size_t blockNumberWithSameName;
	int valueableChainNumber;
	double lengthStandardDeviation;
	double positionStandardDeviation;
	vector<double>* pLengthArray;
	vector<double>* pPositionArray;
	vector<Distribution>* pLengthDistributionArray;
	vector<Distribution>* pPositionDistributionArray;
	BlockFloat() :
		blockName(""),
		blockLength(0),
		blockPosInChain(0),
		blockNumberWithSameName(0),
		valueableChainNumber(0),
		lengthStandardDeviation(0),
		positionStandardDeviation(0),
		pLengthArray(NULL),
		pPositionArray(NULL),
		pLengthDistributionArray(NULL),
		pPositionDistributionArray(NULL)
	{
	}
} BlockFloat;

typedef vector<BlockItem*>* pBlockArray;
typedef vector<string>* pPerChain;
typedef vector<BlockFloat*>* pBlockFloatArray;
typedef vector<Distribution*>* pDistributionArray;
typedef int number_with_sameName;
typedef size_t max_count;

struct fBlockInfo
{
	string _name;
	max_count _maxcount;
	vector<vector<BlockItem*>*>* _pfBlocks;
	fBlockInfo(string name, max_count count, vector<vector<BlockItem*>*>* pV) :_name(name),_maxcount(count),_pfBlocks(pV){}
};

#endif
