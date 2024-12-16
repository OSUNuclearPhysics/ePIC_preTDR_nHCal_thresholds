/*
 * FileList.h
 *
 *  Created on: 1 mar 2022
 *      Author: Khaless
 */

#ifndef FILELIST_H_
#define FILELIST_H_

#include <fstream>
#include <cstring>
#include <string>
#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TChain.h>


using namespace std;

Bool_t addNewFile(TChain *chain, const TString file);
Bool_t openFileList(TChain *chain, const TString listName = "files.list");
std::vector<std::string> openList(TString list);


/////////////////////////////////////////////////////////////
Bool_t addNewFile(TChain *chain, const TString file)
{
	cout << "Opening new data file " << file << " ... ";

	try {
		/* Open the tree file */
		TFile *mFile = TFile::Open(file);
		/* The checks must be in two separate ifs due to underfined processing order otherwise */
		if (mFile == NULL)
			throw kFALSE;
		if (mFile->IsZombie())
			throw kFALSE;
	} catch (...) {
		cout << "FAILURE" << endl;
		return kFALSE;
	}


	cout << "success" << endl;

	chain->Add(file);

	return kTRUE;
}
/////////////////////////////////////////////////////////////
Bool_t openFileList(TChain *chain, const TString listName)
{

	std::ifstream mInputFileList;
	string file;

	mInputFileList.open(listName.Data());
	if (!mInputFileList) {
		cout << "File list " << listName << " couldn't be opened!" << endl;
		return kFALSE;
	}
	cout << "File list " << listName << " open." << endl;


	while(!mInputFileList.eof())	{
		getline(mInputFileList, file);
		if (addNewFile(chain, file) == kFALSE)	continue;
	}

	return kTRUE;
}
/////////////////////////////////////////////////////////////


std::vector<std::string> openList(TString listName)
{
	std::vector<std::string> filenames;

	std::ifstream mInputFileList;
	string file;

	mInputFileList.open(listName.Data());
	if (!mInputFileList) {
		cout << "File list " << listName << " couldn't be opened!" << endl;
		return filenames;
	}
	cout << "File list " << listName << " open." << endl;


	while(!mInputFileList.eof())	{
		getline(mInputFileList, file);
		filenames.push_back(file);
		cout<<"Adding file: "<<file<<endl;

		try {
			/* Open the tree file */
			TFile *mFile = TFile::Open(file.data());
			/* The checks must be in two separate ifs due to underfined processing order otherwise */
			if (mFile == NULL)
				throw kFALSE;
			if (mFile->IsZombie())
				throw kFALSE;
		} catch (...) {
			cout << "FAILURE" << endl;
		}

	}

	return filenames;

}


#endif /* FILELIST_H_ */
