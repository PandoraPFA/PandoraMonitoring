/**
 *  @file   TTreeWrapper/src/TTreeWrapper.cc
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

#include <map>
#include <algorithm>
#include <cctype> // for toupper
#include <cmath>

#include <typeinfo>

namespace pandora_monitoring
{


//------------------------------------------------------------------------------------------------------------------------------------------

TTreeWrapper::TTreeWrapper() 
{
    gROOT->ProcessLine("#include <vector>");
}



//------------------------------------------------------------------------------------------------------------------------------------------

TTreeWrapper::~TTreeWrapper() 
{
    for( TreeMap::iterator itTreeMap = m_treeMap.begin(), itTreeMapEnd = m_treeMap.end(); itTreeMap != itTreeMapEnd; ++itTreeMap )
    {
        delete itTreeMap->second.first; // delete the TTree
        for( BranchMap::iterator itBranchMap = itTreeMap->second.second->begin(), itBranchMapEnd = itTreeMap->second.second->end(); itBranchMap != itBranchMapEnd; ++itBranchMap )
        {
            delete itBranchMap->second; // delete the BranchHandlers
        }
        delete itTreeMap->second.second; // delete the BranchMap itself
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

TTreeWrapper::TreeMap::iterator TTreeWrapper::AddTree( std::string treeName )
{
    TreeMap::iterator treeIt = m_treeMap.find( treeName );
    if( treeIt != m_treeMap.end() )
        return treeIt;

    TTree* tree = new TTree(treeName.c_str(), treeName.c_str() );
    BranchMap* branchMap = new BranchMap();
    std::pair<TreeMap::iterator, bool> itInfo = m_treeMap.insert( TreeMap::value_type( treeName, 
                                                                                       TreeInfo( tree, branchMap ) ) );

    if( !(itInfo.second) )
        throw TreeInsertError();

    treeIt = itInfo.first;
    return treeIt;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template< typename VarType >
TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch( std::string treeName, std::string branchName )
{
    TreeMap::iterator treeIt = m_treeMap.end();
    try
    {
        treeIt = AddTree( treeName );
    }
    catch( TreeInsertError& excpt )
    {
        throw;
    }
    catch(...)
    {
        std::cout << "TTreeWrapper/AddBranch/UNKNOWN EXCEPTION" << std::endl;
        throw;
    }

    BranchMap* branchMap = treeIt->second.second;
    BranchMap::iterator branchIt = branchMap->find( branchName );
    if( branchIt == treeIt->second.second->end() )
    {
        TTree* tree = treeIt->second.first;

        // create a new BranchHandler 
        BranchHandler* branchHandler = new BranchHandler( tree, branchName );

        // insert the branchHandler
        std::pair<BranchMap::iterator,bool> itInfo = treeIt->second.second->insert( BranchMap::value_type( branchName, branchHandler ) );

        if( !itInfo.second )
            throw BranchInsertError();

        return itInfo.first;
    }
    return branchIt;
}


//------------------------------------------------------------------------------------------------------------------------------------------

template< typename VarType >
bool TTreeWrapper::Set( std::string treeName, std::string branchName, VarType value )
{
    BranchMap::iterator branchIt;
    try
    {
        branchIt = AddBranch<VarType>( treeName, branchName );
    }
    catch( BranchInsertError& excpt )
    {
        throw;
    }
    
    BranchHandler* branchHandler = branchIt->second;
    return branchHandler->Set( value );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeWrapper::Fill( std::string treeName )
{
    TreeMap::iterator treeIt = m_treeMap.find( treeName );
    if( treeIt == m_treeMap.end() )
        throw TreeNotFoundError();

    if( treeIt->second.first == NULL )
        throw TreeNotFoundError();

    treeIt->second.first->Fill();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeWrapper::Print( std::string treeName )
{
    TreeMap::iterator treeIt = m_treeMap.find( treeName );
    if( treeIt == m_treeMap.end() )
        throw TreeNotFoundError();

    if( treeIt->second.first == NULL )
        throw TreeNotFoundError();

    treeIt->second.first->Print();
}

//------------------------------------------------------------------------------------------------------------------------------------------

// void TTreeWrapper::Write( std::string treeName )
// {
//     TreeMap::iterator treeIt = m_treeMap.find( treeName );
//     if( treeIt == m_treeMap.end() )
//         throw TreeNotFoundError();

//     if( treeIt->second.first == NULL )
//         throw TreeNotFoundError();

//     treeIt->second.first->Write();
// }

//------------------------------------------------------------------------------------------------------------------------------------------

TTree* TTreeWrapper::GetTree( std::string treeName ) const
{
    TreeMap::const_iterator treeIt = m_treeMap.find( treeName );
    if( treeIt == m_treeMap.end() )
        throw TreeNotFoundError();

    if( treeIt->second.first == NULL )
        throw TreeNotFoundError();

    return treeIt->second.first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TTreeWrapper::Branch<T>::Branch(TTree* tree, std::string branchName) 
    : m_name(branchName), 
      m_pTree(tree), 
      m_pBranch(0) 
{
    const char* name = branchName.c_str();
    std::string typeIdName( typeid(m_variable).name() );
    if( typeIdName.find("vector") != std::string::npos ) // if it is a standard-vector
    {
        isVector = true;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,20,0)
        m_pBranch = tree->Branch( name, &(m_variable)  );
#else
//         m_pBranch = tree->Bronch( name, "float", &(m_variable)  );
#endif
    }
    else
    {
        isVector = false;
        std::transform(typeIdName.begin(), typeIdName.end(), typeIdName.begin(), (int(*)(int))std::toupper);
        if( typeIdName.size() != 1 ) // has to be one letter only to be F, D, I, ...
            throw BadType();
        const char firstLetterOfType = typeIdName.at(0);
        m_pBranch = tree->Branch( name, &m_variable, TString(m_name.c_str())+TString("/")+TString(firstLetterOfType) );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TTreeWrapper::Branch<T>::~Branch()
{
    //    delete m_pBranch;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TTreeWrapper::Branch<T>::Set( T variable )
{
    m_variable = variable;
    if( isVector )
    {
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,20,0)
        m_pBranch->SetAddress( &m_variable );
#endif
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TTreeWrapper::BranchHandler::BranchHandler( TTree* tree, std::string branchName )
    : m_branchType(BRANCH_NO_TYPE_DEFINED), 
      m_branchFloat (NULL),
      m_branchDouble(NULL),
      m_branchInt   (NULL),
      m_branchVectorFloat (NULL),
      m_branchVectorDouble(NULL),
      m_branchVectorInt   (NULL),
      m_tree( tree ),
      m_branchName( branchName )
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

bool TTreeWrapper::BranchHandler::Set( float value )
{
    if( m_branchFloat && m_branchType != BRANCH_FLOAT )
        return false;

    if( m_branchFloat == NULL ) 
    {
        m_branchFloat = new Branch<float>( m_tree, m_branchName ); ///< create a branch of type float
        m_branchType = BRANCH_FLOAT;
    }

    m_branchFloat->Set( value );
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set( double value )
{
    if( m_branchDouble && m_branchType != BRANCH_DOUBLE )
        return false;

    if( m_branchDouble == NULL ) 
    {
        m_branchDouble = new Branch<double>( m_tree, m_branchName ); ///< create a branch of type float
        m_branchType = BRANCH_DOUBLE;
    }

    m_branchDouble->Set( value );
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set( int value )
{
    if( m_branchInt && m_branchType != BRANCH_INT )
        return false;

    if( m_branchInt == NULL ) 
    {
        m_branchInt = new Branch<int>( m_tree, m_branchName ); ///< create a branch of type float
        m_branchType = BRANCH_INT;
    }

    m_branchInt->Set( value );
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set( VectorFloat* ptr )
{
    if( m_branchVectorFloat && m_branchType != BRANCH_VECTOR_FLOAT )
        return false;

    if( m_branchVectorFloat == NULL ) 
    {
        m_branchVectorFloat = new Branch<VectorFloat*>( m_tree, m_branchName ); ///< create a branch of type float
        m_branchType = BRANCH_VECTOR_FLOAT;
    }

    m_branchVectorFloat->Set( ptr );
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set( VectorDouble* ptr )
{
    if( m_branchVectorDouble && m_branchType != BRANCH_VECTOR_DOUBLE )
        return false;

    if( m_branchVectorDouble == NULL ) 
    {
        m_branchVectorDouble = new Branch<VectorDouble*>( m_tree, m_branchName ); ///< create a branch of type Double
        m_branchType = BRANCH_VECTOR_DOUBLE;
    }

    m_branchVectorDouble->Set( ptr );
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TTreeWrapper::BranchHandler::Set( VectorInt* ptr )
{
    if( m_branchVectorInt && m_branchType != BRANCH_VECTOR_INT )
        return false;

    if( m_branchVectorInt == NULL ) 
    {
        m_branchVectorInt = new Branch<VectorInt*>( m_tree, m_branchName ); ///< create a branch of type int
        m_branchType = BRANCH_VECTOR_INT;
    }

    m_branchVectorInt->Set( ptr );
    return true;
}


// member function template initializations
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<float>( std::string, std::string );
template bool TTreeWrapper::Set<float>( std::string, std::string, float );
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<int>( std::string, std::string );
template bool TTreeWrapper::Set<int>( std::string, std::string, int );
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<double>( std::string, std::string );
template bool TTreeWrapper::Set<double>( std::string, std::string, double );

template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<VectorFloat*>( std::string, std::string );
template bool TTreeWrapper::Set<VectorFloat*>( std::string, std::string, VectorFloat* );
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<VectorInt*>( std::string, std::string );
template bool TTreeWrapper::Set<VectorInt*>( std::string, std::string, VectorInt* );
template TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch<VectorDouble*>( std::string, std::string );
template bool TTreeWrapper::Set<VectorDouble*>( std::string, std::string, VectorDouble* );




} // namespace pandora_monitoring
