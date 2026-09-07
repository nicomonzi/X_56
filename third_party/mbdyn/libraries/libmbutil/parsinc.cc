/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* parser con include */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstring>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_PWD_H
#include <pwd.h>
#endif /* HAVE_PWD_H */
#include <errno.h>

#include "parsinc.h"
#include "filename.h"

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif // !PATH_MAX



struct IncludeDR : public DescRead {
	bool Read(HighParser& HP);
};

bool
IncludeDR::Read(HighParser& HP)
{
	IncludeParser *pIP = dynamic_cast<IncludeParser *>(&HP);
	ASSERT(pIP != 0);
	return pIP->Include_int();
}

struct ChDirDR : public DescRead {
	bool Read(HighParser& HP);
};

bool
ChDirDR::Read(HighParser& HP)
{
	if (!HP.IsArg()) {
		silent_cerr("Parser error "
			"in ChDir::Read(), "
			"colon expected at line " << HP.GetLineData() 
			<< std::endl);
      		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
   	}

	IncludeParser *pIP = dynamic_cast<IncludeParser *>(&HP);
	ASSERT(pIP != 0);
   
	const std::string sfname = pIP->GetFileName();

	if (chdir(sfname.c_str())) {
		silent_cerr("Error in chdir, path=\"" << sfname << "\" at line "
			<< HP.GetLineData() << std::endl);
		throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
	}

	return true;
}

static unsigned desc_done;

static void
InitDescData(void)
{
	// NOTE: data will be destroyed when the underlying HighParser is destroyed (is this what we want?)
	if (::desc_done++ > 0) {
		return;
	}

	SetDescData("include", new IncludeDR);
	SetDescData("chdir", new ChDirDR);

	/* NOTE: add here initialization of new built-in descriptions;
	 * alternative ways to register new custom descriptions are:
	 * - call SetDescData() from anywhere in the code
	 * - write a module that calls SetDescData() from inside a function
	 *   called module_init(), and run-time load it using "module load"
	 *   in the input file.
	 */
}


/* IncludeParser - begin */

IncludeParser::IncludeParser(MathParser& MP,
			     InputStream& streamIn,
			     const std::string sInitialFile)
: HighParser(MP, streamIn)
{
	ASSERT(!sInitialFile.empty());
#ifdef USE_INCLUDE_PARSER
   	char s[PATH_MAX];
   	if (getcwd(s, sizeof(s)) == NULL) {
		silent_cerr("Error in getcwd()" << std::endl);
      		throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
   	}
	sCurrPath = s;
	sInitialPath = sCurrPath;
   	DEBUGCOUT("Current directory is \"" << sCurrPath << "\"" << std::endl);

   	sCurrFile = sInitialFile;
#else /* !USE_INCLUDE_PARSER */
   	NO_OP;
#endif /* !USE_INCLUDE_PARSER */

	// NOTE: data will be destroyed when the underlying HighParser is destroyed (is this what we want?)
	InitDescData();
}


IncludeParser::~IncludeParser(void)
{   
   	IncludeParser::Close();
}
 

void IncludeParser::Close(void)
{
   	MyInput* pmi = NULL;
	if (!myinput.empty()) {
		pmi = myinput.top();
      		ASSERT(pmi != NULL);
      		/* Nota: deve esserci solo l'ultimo file */
      		ASSERT(pf != NULL);
      		ASSERT(pIn != NULL);

#ifdef USE_INCLUDE_PARSER
      		ASSERT(!sCurrPath.empty());
      		ASSERT(!sCurrFile.empty());
#endif /* USE_INCLUDE_PARSER */

      		if (pf != NULL) {
	 		SAFEDELETE(pf);
      		}
      		if (pIn != NULL) {
	 		SAFEDELETE(pIn);
      		}

#ifdef USE_INCLUDE_PARSER
      		DEBUGCOUT("Leaving directory <" << sCurrPath
			<< ">, file <" << sCurrFile << '>' << std::endl);
#endif /* USE_INCLUDE_PARSER */

      		pf = pmi->pfile;
      		pIn = pmi->pis;

#ifdef USE_INCLUDE_PARSER
      		sCurrPath = pmi->sPath;
      		sCurrFile = pmi->sFile;
      		DEBUGCOUT("Entering directory \"" << sCurrPath
			<< "\", file \"" << sCurrFile << "\"" << std::endl);
      		if (chdir(sCurrPath.c_str())) {
			silent_cerr("Error in chdir, path=\""
				<< sCurrPath << "\"" << std::endl);
	 		throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
      		}
#endif /* USE_INCLUDE_PARSER */

      		/* pmi must be non NULL */
      		SAFEDELETE(pmi);

		myinput.pop();
   	}

#ifdef USE_INCLUDE_PARSER
   	sCurrPath.clear();
   	sCurrFile.clear();
#endif /* USE_INCLUDE_PARSER */
}

flag 
IncludeParser::fCheckStack(void)
{
   	MyInput* pmi = NULL;
	if (!myinput.empty()) {
		pmi = myinput.top();

      		ASSERT(pmi != NULL);
      		/* 
       		 * Nota: se la stack e' piena, allora sia pf che pIn
		 * devono essere diversi da NULL; viceversa, se la stack
		 * e' vuota, pf deve essere NULL.
		 */
      		ASSERT(pf != NULL);
      		ASSERT(pIn != NULL);
#ifdef USE_INCLUDE_PARSER
      		ASSERT(!sCurrPath.empty());
      		ASSERT(!sCurrFile.empty());
#endif /* USE_INCLUDE_PARSER */

      		SAFEDELETE(pf);
      		SAFEDELETE(pIn);
#ifdef USE_INCLUDE_PARSER
      		DEBUGCOUT("Leaving directory <" << sCurrPath
			<< ">, file <" << sCurrFile << '>' << std::endl);
#endif /* USE_INCLUDE_PARSER */

      		pf = pmi->pfile;
      		pIn = pmi->pis;
#ifdef USE_INCLUDE_PARSER
      		sCurrPath = pmi->sPath;
      		sCurrFile = pmi->sFile;
      		DEBUGCOUT("Entering directory \"" << sCurrPath
			<< "\", file \"" << sCurrFile << "\"" << std::endl);
      		if (chdir(sCurrPath.c_str())) {
			silent_cerr("Error in chdir, path=\""
				<< sCurrPath << "\"" << std::endl);
	 		throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
      		}
#endif /* USE_INCLUDE_PARSER */
      
      		SAFEDELETE(pmi);

		myinput.pop();

      		return flag(1);

   	} else {
      		return flag(0);
   	}
}

bool
IncludeParser::Include_int()
{
   	// if (FirstToken() == UNKNOWN) {
   	if (!IsArg()) {
		silent_cerr("Parser error in IncludeParser::Include_int(),"
			" colon expected at line " << GetLineData() 
			<< std::endl);
      		throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
   	}
   
	const std::string sfname = GetFileName();

	if (!sfname.empty()) {
		struct stat	s;

		if (stat(sfname.c_str(), &s)) {
			int save_errno = errno;

			silent_cerr("Cannot stat file <" << sfname << "> "
				"at line " << GetLineData() << ": "
				<< save_errno 
				<< " (" << strerror(save_errno) << ")" 
				<< std::endl);
			throw ErrFile(MBDYN_EXCEPT_ARGS);
		}

		if (!S_ISREG(s.st_mode)) {
			silent_cerr("File <" << sfname << "> "
				"at line " << GetLineData() << ": "
				"not a regular file?" << std::endl);
			throw ErrFile(MBDYN_EXCEPT_ARGS);
		}

		if (!(s.st_mode & S_IRUSR)) {
			silent_cerr("File <" << sfname << "> "
				"at line " << GetLineData() << ": "
				"no read permissions?" << std::endl);
			throw ErrFile(MBDYN_EXCEPT_ARGS);
		}
	} else {
		silent_cerr("File name expected at line " << GetLineData() << std::endl);
		throw ErrFile(MBDYN_EXCEPT_ARGS);
	}

	std::ifstream *pf_old = pf;
	InputStream *pIn_old = pIn;
	std::string sOldPath = sCurrPath;
	std::string sOldFile = sCurrFile;

   	pf = NULL;
   	pIn = NULL;

#ifdef _WIN32
	// open the file in non translated mode in order not to break seek operations
	SAFENEWWITHCONSTRUCTOR(pf, std::ifstream, std::ifstream(sfname, std::ios::binary));
#else
	SAFENEWWITHCONSTRUCTOR(pf, std::ifstream, std::ifstream(sfname));
#endif
   	if (!(*pf)) {
#ifdef DEBUG
		char *buf = getcwd(NULL, 0);
		if (buf != NULL) {
			DEBUGCERR("Current directory \"" << buf << "\"" 
					<< std::endl);
			free(buf);
		}
#endif /* DEBUG */

		/* restore */
		pf = pf_old;
		pIn = pIn_old;
		sCurrPath = sOldPath;
		sCurrFile = sOldFile;
   
		silent_cerr("Invalid file <" << sfname << "> "
			"at line " << GetLineData() << std::endl);
      		throw ErrFile(MBDYN_EXCEPT_ARGS);
   	}
   
   	SAFENEWWITHCONSTRUCTOR(pIn, InputStream, InputStream(*pf));

   	/* Cambio di directory */
#ifdef USE_INCLUDE_PARSER
   	sCurrPath.clear();
   	sCurrFile.clear();

	std::string::size_type sep = sfname.find_last_of(DIR_SEP);
	if (sep != std::string::npos) {
		std::string sDir(sfname, 0, sep + 1);
 		if (chdir(sDir.c_str())) {
			silent_cerr("Error in chdir, path="
				<< sDir << std::endl);
    			throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
 		}
 		char p[PATH_MAX];
 		if (getcwd(p, sizeof(p)) == NULL) {
			silent_cerr("Error in getcwd()" << std::endl);
    			throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
 		}
		sCurrPath = p;
 		DEBUGCOUT("Current directory is \"" << sCurrPath
			<< "\"" << std::endl);

		sCurrFile = sfname.substr(sep + 1);

	} else {
		sCurrFile = sfname;
	}
   	DEBUGCOUT("Opening file <" << sCurrFile << '>' << std::endl);

   	if (sCurrPath.empty()) {
      		char ss[PATH_MAX];
      		if (getcwd(ss, sizeof(ss)) == NULL) {
			silent_cerr("Error in getcwd()" << std::endl);
	 		throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
      		}
		sCurrPath = ss;
      		DEBUGCOUT("Current directory is \"" << sCurrPath
			<< "\"" << std::endl);
   	}
#endif /* USE_INCLUDE_PARSER */

      	MyInput* pmi = NULL;
#ifdef USE_INCLUDE_PARSER
   	SAFENEWWITHCONSTRUCTOR(pmi, 
			       MyInput,
			       MyInput(pf_old, pIn_old, sOldPath, sOldFile));
#else /* !USE_INCLUDE_PARSER */
   	SAFENEWWITHCONSTRUCTOR(pmi, MyInput, MyInput(pf_old, pIn_old));
#endif /* !USE_INCLUDE_PARSER */
 
	myinput.push(pmi);
 
   	/* FIXME: mettere un test se c'e' il punto e virgola? */
   	CurrToken = HighParser::DESCRIPTION;

	return true;
}

void
IncludeParser::Eof(void)
{
	if (!fCheckStack()) {
		throw EndOfFile(MBDYN_EXCEPT_ARGS);
	}
}

/*
 * expands environment variables in "in" into "out";
 * returns false (leaving "out" in an undefined state) on failure
 */
static bool
expand_environment(const std::string& in, std::string& out)
{
	DEBUGCOUT(">> expand_environment: " << in << std::endl);

	out.clear();
	for (std::string::size_type c = 0; c < in.size(); c++) {
		if (in[c] != '$') {
			out += in[c];
			continue;
		}

		/* "$$" is a literal '$' */
		if (in[c + 1] == '$') {
			out += '$';
			c++;
			continue;
		}

		c++;
		std::string::size_type namepos = c;
		const char *value = NULL;
		if (in[c] == '{') {
			const std::string::size_type end = in.find('}', c);

			if (end == std::string::npos) {
				silent_cerr("missing trailing \"}\" "
						"in \"" << in << "\""
						<< std::endl);
				return false;
			}

			namepos++;
			std::string buf(in, namepos, end - namepos);
			value = getenv(buf.c_str());
			if (value == NULL) {
				silent_cerr("unable to find "
						"environment "
						"variable \""
						<< buf << "\""
						<< std::endl);
				return false;
			}

			/* skip past the closing brace ('}';
			 * the for loop increments c) */
			c = end;

		} else {
			if (in[c] != '_' && !isalpha(in[c])) {
				silent_cerr("illegal leading char "
						"in environment "
						"variable name in \""
						<< in << "\""
						<< std::endl);
				return false;
			}

			for (c++; c < in.size(); c++) {
				if (in[c] != '_' && !isalnum(in[c])) {
					break;
				}
			}

			std::string buf(in, namepos, c - namepos);
			value = getenv(buf.c_str());
			if (value == NULL) {
				silent_cerr("unable to find "
						"environment "
						"variable \""
						<< buf << "\""
						<< std::endl);
				return false;
			}

			/* because it's incremented again by "for" */
			c--;
		}

		out += value;
	}

	DEBUGCOUT("<< expand_environment: " << out << std::endl);

	return true;
}

/*
 * resolves environment variables and "~" prefixes in "filename_in";
 * returns false (leaving "res" in an undefined state) on failure
 */
static bool
resolve_filename(const std::string& filename_in, std::string& res)
{
        std::string filename;

        if (filename_in.find('$') != std::string::npos) {
                if (!expand_environment(filename_in, filename)) {
                        return false;
                }

        } else {
                filename = filename_in;
        }

        if (!filename.empty() && filename[0] == '~') {
                if (filename.length() > 1 && filename[1] == DIR_SEP) {
                        /* do environment stuff */
                        const char *home = getenv("HOME");
                        if (home != NULL) {
                                res = home + filename.substr(1);
                                return true;
                        }

#if defined(HAVE_PWD_H)
                } else {
                        std::string::size_type p = filename.find(DIR_SEP, 1);
                        if (p != std::string::npos) {
                                std::string buf = filename.substr(1, p - 1);

                                /* do passwd stuff */
                                struct passwd *pw = getpwnam(buf.c_str());
                                if (pw != NULL) {
                                        res = pw->pw_dir + filename.substr(p);
                                        return true;
                                }
                        }
#endif /* HAVE_PWD_H */
                }
        }

        res = filename;

        return true;
}

std::string
IncludeParser::GetFileName(enum Delims Del)
{
	const std::string s = GetStringWithDelims(Del);
	if (s.empty()) {
		return std::string();
	}

   	std::string stmp;
	if (!resolve_filename(s, stmp)) {
		return std::string();
	}

	return stmp;
}

HighParser::ErrOut
IncludeParser::GetLineData(void) const
{
   	ErrOut LineData;
	LineData.sFileName = sCurrFile;
	LineData.sPathName = (sCurrPath.empty() || sCurrPath == sInitialPath)
		? std::string() : sCurrPath;
   	LineData.iLineNumber = GetLineNumber();
   	return LineData;
}

/* IncludeParser - end */
