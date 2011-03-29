/**
 *  @file   src/TTreeWrapper.cc
 * 
 *  @brief  Implementation of the pandora monitoring class.
 * 
 *  $Log: $
 */

// ROOT include files
#include "TSystem.h"
#include "TTree.h"
#include "TROOT.h" 
#include "TBranch.h"

#include "TTreeWrapper.h"

#include <algorithm>
#include <cctype> // for toupper
#include <cmath>
#include <typeinfo>

namespace pandora_monitoring
{

TTreeWrapper::TTreeWrapper() 
{
    gROOT->ProcessLine("#include <vector>"); // should not be necessary for newer versions of ROOT
}

//------------------------------------------------------------------------------------------------------------------------------------------

TTreeWrapper::~TTreeWrapper() 
{
    Clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeWrapper::Clear() 
{
    for(TreeMap::iterator itTreeMap = m_treeMap.begin(), itTreeMapEnd = m_treeMap.end(); itTreeMap != itTreeMapEnd; ++itTreeMap)
    {
        TTree *tree = itTreeMap->second.first;
        if (tree)
        {
            delete itTreeMap->second.first;     // delete the TTree
            itTreeMap->second.first = NULL;

            for(BranchMap::iterator itBranchMap = itTreeMap->second.second->begin(), itBranchMapEnd = itTreeMap->second.second->end();
                itBranchMap != itBranchMapEnd; ++itBranchMap)
            {
                delete itBranchMap->second;     // delete the BranchHandlers
            }
        }

        delete itTreeMap->second.second;        // delete the BranchMap itself
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template< typename VarType >
bool TTreeWrapper::Set(const std::string &treeName, const std::string &branchName, VarType value)
{
    BranchMap::iterator branchIt;

    try
    {
        branchIt = AddBranch<VarType>(treeName, branchName);
    }
    catch(BranchInsertError &excpt)
    {
        throw;
    }

    BranchHandler *branchHandler = branchIt->second;

    return branchHandler->Set(value);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeWrapper::Fill(const std::string &treeName)
{
    TreeMap::iterator treeIt = m_treeMap.find(treeName);

    if(treeIt == m_treeMap.end())
        throw TreeNotFoundError();

    if(treeIt->second.first == NULL)
        throw TreeNotFoundError();

    treeIt->second.first->Fill();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeWrapper::Print(const std::string &treeName) const
{
    TreeMap::const_iterator treeIt = m_treeMap.find(treeName);

    if(treeIt == m_treeMap.end())
        throw TreeNotFoundError();

    if(treeIt->second.first == NULL)
        throw TreeNotFoundError();

    treeIt->second.first->Print();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeWrapper::Scan(const std::string &treeName) const
{
    TreeMap::const_iterator treeIt = m_treeMap.find(treeName);

    if(treeIt == m_treeMap.end())
        throw TreeNotFoundError();

    if(treeIt->second.first == NULL)
        throw TreeNotFoundError();

     bool isFirst = true;
     std::string branchNames;

     for (BranchMap::iterator itBranchMap = treeIt->second.second->begin(), itBranchMapEnd = treeIt->second.second->begin();
        itBranchMap != itBranchMapEnd; ++itBranchMap)
     {
         if(!isFirst)
         {
             branchNames += ":";
         }

         branchNames += itBranchMap->first;
         isFirst = false;
     }

     treeIt->second.first->Scan(branchNames.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

TTree *&TTreeWrapper::GetTree(const std::string &treeName)
{
    TreeMap::iterator treeIt = m_treeMap.find(treeName);

    if(treeIt == m_treeMap.end())
        throw TreeNotFoundError();

    if(treeIt->second.first == NULL)
        throw TreeNotFoundError();

    return treeIt->second.first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TTreeWrapper::TreeMap::iterator TTreeWrapper::AddTree(const std::string &treeName)
{
    TreeMap::iterator treeIt = m_treeMap.find(treeName);

    if(treeIt != m_treeMap.end())
        return treeIt;

    TTree* tree = new TTree(treeName.c_str(), treeName.c_str());

    BranchMap* branchMap = new BranchMap();
    std::pair<TreeMap::iterator, bool> itInfo = m_treeMap.insert(TreeMap::value_type(treeName, TreeInfo(tree, branchMap)));

    if(!(itInfo.second))
        throw TreeInsertError();

    treeIt = itInfo.first;

    return treeIt;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template< typename VarType >
TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch(const std::string &treeName, const std::string &branchName)
{
    TreeMap::iterator treeIt = m_treeMap.end();

    try
    {
        treeIt = AddTree(treeName);
    }
    catch (TreeInsertError& excpt)
    {
        throw;
    }
    catch (...)
    {
        std::cout << "TTreeWrapper/AddBranch/UNKNOWN EXCEPTION" << std::endl;
        throw;
    }

    BranchMap::iterator branchIt = treeIt->second.second->find(branchName);

    if (branchIt == treeIt->second.second->end())
    {
        TTree *tree = treeIt->second.first;

        // create a new BranchHandler 
        BranchHandler* branchHandler = new BranchHandler(tree, branchName);

        // insert the branchHandler
        std::pair<BranchMap::iterator,bool> itInfo = treeIt->second.second->insert(BranchMap::value_type(branchName, branchHandler));

        if (!itInfo.second)
            throw BranchInsertError();

        return itInfo.first;
    }

    return branchIt;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TTreeWrapper::Branch<T>::Branch(TTree *pTree, const std::string &branchName) :
    m_name(branchName), 
    m_pTree(pTree), 
    m_pBranch(0) 
{
    const char* name = branchName.c_str();
    std::string typeIdName(typeid(m_variable).name());

    if (typeIdName.find("vector") != std::string::npos)
    {
        m_isVector = true;

        m_pBranch = pTree->Branch(name, &m_variable );

        if (!m_pBranch)
        {
            std::cout << "TTreeWrapper::Branch<T>::Branch(...) with [ T = " << typeid( T ).name() << " ], m_pBranch = " << m_pBranch
                      << " : branch has not been created." << std::endl;
            throw BranchInsertError();
        }
    }
    else
    {
        m_isVector = false;
        std::transform(typeIdName.begin(), typeIdName.end(), typeIdName.begin(), (int(*)(int))std::toupper);

        if (typeIdName.size() != 1) // has to be one letter only to be F, D, I, ...
            throw BadType();

        const char firstLetterOfType = typeIdName.at(0);
        m_pBranch = pTree->Branch(name, &m_variable, TString(m_name.c_str()) + TString("/") + TString(firstLetterOfType));

        if (!m_pBranch)
        {
            std::cout << "TTreeWrapper::Branch<T>::Branch(...) with [ T = " << typeid( T ).name() << " ], m_pBranch = " << m_pBranch
                      << " : branch has not been created." << std::endl;
            throw BranchInsertError();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TTreeWrapper::Branch<T>::~Branch()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TTreeWrapper::Branch<T>::Set(T variable)
{
    m_variable = variable;

    if(m_isVector)
    {
        if( !m_pBranch )
        {
            std::cout << "TTreeWrapper::Branch<T>::Set( T variable ) with [ T = " << typeid( T ).name() << " ], m_pBranch = " << m_pBranch << std::endl;
            throw BranchNotFoundError();
        }
        m_pBranch->SetAddress(&m_variable);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TTreeWrapper::BranchHandler::BranchHandler(TTree *pTree, const std::string &branchName) :
    m_branchType(BRANCH_NO_TYPE_DEFINED), 
    m_branchFloat (NULL),
    m_branchDouble(NULL),
    m_branchInt   (NULL),
    m_branchVectorFloat (NULL),
    m_branchVectorDouble(NULL),
    m_branchVectorInt   (NULL),
    m_tree(pTree),
    m_branchName(branchName)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TTreeWrapper::BranchHandler::~BranchHandler()
{
    delete m_branchFloat;
    delete m_branchDouble;
    delete m_branchInt;
    delete m_branchVectorFloat;
    delete m_branchVectorDouble;
    delete m_branchVectorInt;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set(float value)
{
    if(m_branchFloat && m_branchType != BRANCH_FLOAT)
        return false;

    if(m_branchFloat == NULL) 
    {
        m_branchFloat = new Branch<float>(m_tree, m_branchName); ///< create a branch of type float
        m_branchType = BRANCH_FLOAT;
    }

    m_branchFloat->Set(value);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set(double value)
{
    if (m_branchDouble && m_branchType != BRANCH_DOUBLE)
        return false;

    if(m_branchDouble == NULL) 
    {
        m_branchDouble = new Branch<double>(m_tree, m_branchName); ///< create a branch of type double
        m_branchType = BRANCH_DOUBLE;
    }

    m_branchDouble->Set(value);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set(int value)
{
    if (m_branchInt && m_branchType != BRANCH_INT)
        return false;

    if (m_branchInt == NULL) 
    {
        m_branchInt = new Branch<int>(m_tree, m_branchName); ///< create a branch of type int
        m_branchType = BRANCH_INT;
    }

    m_branchInt->Set(value);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set(VectorFloat *ptr)
{
    if(m_branchVectorFloat && m_branchType != BRANCH_VECTOR_FLOAT)
        return false;

    if(m_branchVectorFloat == NULL) 
    {
        m_branchVectorFloat = new Branch<VectorFloat*>(m_tree, m_branchName); ///< create a branch of type VectorFloat
        m_branchType = BRANCH_VECTOR_FLOAT;
    }

    m_branchVectorFloat->Set(ptr);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set(VectorDouble *ptr)
{
    if(m_branchVectorDouble && m_branchType != BRANCH_VECTOR_DOUBLE)
        return false;

    if(m_branchVectorDouble == NULL) 
    {
        m_branchVectorDouble = new Branch<VectorDouble*>(m_tree, m_branchName); ///< create a branch of type VectorDouble
        m_branchType = BRANCH_VECTOR_DOUBLE;
    }

    m_branchVectorDouble->Set(ptr);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set(VectorInt *ptr)
{
    if(m_branchVectorInt && m_branchType != BRANCH_VECTOR_INT)
        return false;

    if(m_branchVectorInt == NULL) 
    {
        m_branchVectorInt = new Branch<VectorInt*>(m_tree, m_branchName); ///< create a branch of type VectorInt

        m_branchType = BRANCH_VECTOR_INT;
    }

    m_branchVectorInt->Set(ptr);

    return true;
}

// member function template initializations
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<float>(const std::string &, const std::string &);
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<int>(const std::string &, const std::string &);
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<double>(const std::string &, const std::string &);
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<VectorFloat*>(const std::string &, const std::string &);
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<VectorDouble*>(const std::string &, const std::string &);
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<VectorInt*>(const std::string &, const std::string &);

template bool TTreeWrapper::Set<float>(const std::string &, const std::string &, float);
template bool TTreeWrapper::Set<int>(const std::string &, const std::string &, int);
template bool TTreeWrapper::Set<double>(const std::string &, const std::string &, double);
template bool TTreeWrapper::Set<VectorFloat*>(const std::string &, const std::string  &, VectorFloat*);
template bool TTreeWrapper::Set<VectorInt*>(const std::string &, const std::string &, VectorInt*);
template bool TTreeWrapper::Set<VectorDouble*>(const std::string &, const std::string &, VectorDouble*);

} // namespace pandora_monitoring
