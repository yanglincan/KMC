#include "reaction.h"
using namespace std;

Reaction::~Reaction() 
{
	DeleteElementVector();
}

ostream& operator << (ostream& os, const Reaction& reaction)
{
	if (!reaction.IsReactionSucceeded())
	{
		os << "No reaction!" << endl;
		return os;
	}
	for (auto itor = reaction.m_pElementVector->begin(); itor != reaction.m_pElementVector->end(); itor++)
	{
		os << (*itor)->stoichiometricNumber << "[" << (*itor)->elementName << "]";
		REACTANT == (*itor)->elementType ? os << "{Reactant, " : os << "{Product, ";
		switch ((*itor)->elementRole)
		{
		case UNSPECIFIED: {os << "Unspecified}"; break; }
		case MONOMER: {os << "Monomer}"; break; }
		case SPECIE: {os << "Specie}"; break; }
		case INITIATOR: {os << "Initiator}"; break; }
		}
		os << "  ";
	}
	os << "  rate: " << reaction.m_ReactionRateConstant << 
		"  order("<<reaction.m_RelationBetween_kMC_k.x<<','<< reaction.m_RelationBetween_kMC_k.y<<")"<<endl;
	return os;
}

void Reaction::DeleteElementVector()
{
	if (NULL == m_pElementVector)
		return;
	for (auto itor = m_pElementVector->begin(); itor != m_pElementVector->end(); itor++)
		delete *itor;
	delete m_pElementVector;
	m_pElementVector = NULL;
}

Reaction::Reaction(string strReaction):
	m_pElementVector(NULL), 
	m_ReactionOrder(_UNSPECIFIED_TYPE),
	m_ReactionRateConstant(0),
	m_RelationBetween_kMC_k(OrderController(0,0))
{
	try {
		CheckSyntaxAndInitReac(strReaction);
		CheckReactionOrder();
	}
	catch (exception & e){
		cout << e.what() << endl;
	}
}

Reaction& Reaction::operator = (const Reaction& reaction)
{
	if (!reaction.IsReactionSucceeded())
		return *this;
	this->m_ReactionOrder = reaction.GetReactionOrder();
	this->m_RelationBetween_kMC_k = reaction.GetRelationBetween_kMC_k();
	this->m_ReactionRateConstant = reaction.GetReactionRateConstant();
	this->m_pElementVector = new vector<ReactionElement*>;
	for (auto itor = reaction.GetElementVector().begin(); itor != reaction.GetElementVector().end(); itor++)
		m_pElementVector->push_back(new ReactionElement(
		(*itor)->elementName,
			(*itor)->elementType,
			(*itor)->elementRole,
			(*itor)->stoichiometricNumber));
	return *this;
}

bool BracketsIsPair(char left, char right)
{
	return ('(' == left && ')' == right || '[' == left && ']' == right || '{' == left && '}' == right) ? 
		true : false;
}

//括号配对检测  ()、[]、{}
bool BracketCheck(string str, int low, int high)
{
	stack<char> S;	//使用栈记录已发现但尚未匹配的左括号
	for (int i = low; i < high; i++)
		if ('(' == str[i] || '[' == str[i] || '{' == str[i])
			S.push(str[i]);
		else if (')' == str[i] || ']' == str[i] || '}' == str[i])
				if (!S.empty())
					if (true == BracketsIsPair(S.top(), str[i])) S.pop();
					else 
						return false;
				else 
					return false;
	return S.empty();
}

//删除s字符串中所有的mark子串
void DeleteAllMark(string& s, const string& mark)
{
	size_t nSize = mark.size();
	while (1)
	{
		size_t pos = s.find(mark);
		if (pos == string::npos)
			return;
		s.erase(pos, nSize);
	}
}

void Reaction::CheckSyntaxAndInitReac(string strReaction)
{
	if (strReaction.empty())
		ThrowExceptionAndDestruction(NO_REACTION);
	DeleteAllMark(strReaction, " ");	//删除所有空格
	if(!BracketCheck(strReaction, 0, (int)strReaction.size()))
		ThrowExceptionAndDestruction(BRACKET_ERROR);
	if (string::npos == strReaction.find('('))
		ThrowExceptionAndDestruction(NON_SMALL_BRACKET);
	if (string::npos == strReaction.find('['))
		ThrowExceptionAndDestruction(NON_BRACKET);
	if (string::npos == strReaction.find('>'))
		ThrowExceptionAndDestruction(NON_ARROW);
	if (string::npos == strReaction.find('{'))
		ThrowExceptionAndDestruction(NON_BIG_BRACKET);
	vector<string> strVector;
	size_t size = strReaction.size();
	size_t leftPos = 0;
	for (size_t i = 0; i < size; i++)
	{
		if ('(' == strReaction[i] || '[' == strReaction[i] || '{' == strReaction[i])
			leftPos = i;
		else if (')' == strReaction[i] || ']' == strReaction[i] || '}' == strReaction[i])
			strVector.push_back(strReaction.substr(leftPos + 1, i - leftPos - 1));
		else if ('>' == strReaction[i])
			strVector.push_back(string(1,strReaction[i]));
	}
	const size_t sizeVec = strVector.size();
	bool arrowFlag(false), reactantFlag(false), productFlag(false);
	m_pElementVector = new vector<ReactionElement*>;
	for(size_t i = 0; i < sizeVec; i++)
	{
		if (">" == strVector[i])
		{
			arrowFlag = true;
			continue;
		}
		if (sizeVec-1 == i)
		{
			m_ReactionRateConstant = atof(strVector[i].data());
			if (0 == m_ReactionRateConstant || 0 > m_ReactionRateConstant)
				ThrowExceptionAndDestruction(NON_RATE);
		}
		else
		{
			SubtanceType substanceType = arrowFlag ? PRODUCT : REACTANT;
			REACTANT == substanceType ? reactantFlag = true : productFlag = true;
			int number = atoi(strVector[i].data());
			if (0 == number || 0 > number)
				ThrowExceptionAndDestruction(NON_STOICHIOMETRIC_NUMBER);
			string name = strVector[i + 1];
			if(name.empty())
				ThrowExceptionAndDestruction(NON_NAME);
			m_pElementVector->push_back(new ReactionElement(name, substanceType, UNSPECIFIED, number));
			i++;
		}
	}
	if(!reactantFlag || !productFlag)
		ThrowExceptionAndDestruction(ONLY_REACTANT_OR_PRODUCT);
}

void Reaction::ZeroReactant()
{	
	m_ReactionOrder = _ZERO_TYPE;
	m_RelationBetween_kMC_k = OrderController(-1, -1);
}
void Reaction::OneReactant()
{
	m_ReactionOrder = _A_TYPE;
	m_RelationBetween_kMC_k = OrderController(1, 0);
}
void Reaction::TwoReactant()
{
	ReactionElement* pReactant[2] = {NULL};
	int i = 0;
	for (auto itor = m_pElementVector->begin(); itor != m_pElementVector->end(); itor++)
	{
		if (REACTANT == (*itor)->elementType)
			pReactant[i++] = *itor;
	}
	if ((NULL != pReactant[0] && NULL == pReactant[1]) ||
		(pReactant[0]->elementName == pReactant[1]->elementName))
	{
		m_ReactionOrder = _2A_TYPE;
		m_RelationBetween_kMC_k = OrderController(2, 1);
	}
	else
	{
		m_ReactionOrder = _AB_TYPE;
		m_RelationBetween_kMC_k = OrderController(1, 1);
	}
}
void Reaction::ThreeReactant()
{
	ReactionElement* pReactant[3] = { NULL };
	int i = 0;
	for (auto itor = m_pElementVector->begin(); itor != m_pElementVector->end(); itor++)
	{
		if (REACTANT == (*itor)->elementType)
			pReactant[i++] = *itor;
	}
	if ((NULL != pReactant[0] && NULL == pReactant[1] && NULL == pReactant[2]) ||
		(pReactant[0]->elementName == pReactant[1]->elementName && pReactant[1]->elementName == pReactant[2]->elementName))
	{
		m_ReactionOrder = _3A_TYPE;
		m_RelationBetween_kMC_k = OrderController(6, 2);
	}
	else if ((NULL != pReactant[0] && NULL != pReactant[1] && NULL == pReactant[2]) ||
			(pReactant[0]->elementName == pReactant[1]->elementName && pReactant[1]->elementName != pReactant[2]->elementName) ||
			(pReactant[1]->elementName == pReactant[2]->elementName && pReactant[2]->elementName != pReactant[0]->elementName) ||
			(pReactant[2]->elementName == pReactant[0]->elementName && pReactant[0]->elementName != pReactant[1]->elementName))
	{
		m_ReactionOrder = _A2B_YTPE;
		m_RelationBetween_kMC_k = OrderController(2, 2);
	}
	else
	{
		m_ReactionOrder = _ABC_TYPE;
		m_RelationBetween_kMC_k = OrderController(1, 2);
	}
}

void Reaction::MoreThanThreeReactant()
{
	ThrowExceptionAndDestruction(MORE_THAN_THREE);
}

void Reaction::CheckReactionOrder()
{
	if (_UNSPECIFIED_TYPE != m_ReactionOrder)
		return;

	int nNum_Reactant = 0;
	for (auto itor = m_pElementVector->begin(); itor!= m_pElementVector->end(); itor++)
		if (REACTANT == (*itor)->elementType)
			nNum_Reactant += (*itor)->stoichiometricNumber;

	switch (nNum_Reactant)
	{
		case 0: ZeroReactant(); break;
		case 1: OneReactant(); break; 
		case 2: TwoReactant(); break; 
		case 3: ThreeReactant(); break;
		default: MoreThanThreeReactant();
	}
}