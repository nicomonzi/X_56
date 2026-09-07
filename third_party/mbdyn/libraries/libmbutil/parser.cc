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

/* parser */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstring>
#include <stdlib.h>
#include <stack>
#include <map>
#include <list>

#include "mathtyp.h"
#include "parser.h"
#include "filename.h"
#include "Rot.hh"

/* LowParser - begin */

static int
skip_remarks(HighParser& HP, InputStream& In, char &cIn)
{
skip_again:;
	for (;;) {
		In.get(cIn);
                if (In.eof()) {
                        return -1;
                }
		if (!isspace(cIn)) {
			break;
		}
	}

        switch (cIn) {
        case MathParser::ONE_LINE_REMARK:
		for (;;) {
			In.get(cIn);
			if (In.eof()) {
				return -1;
			}
			if (cIn == '\n') {
				break;
			}
                        if (cIn == '\\') {
                                In.get(cIn);
                                if (In.eof()) {
                                        return -1;
                                }
                                if (cIn == '\r') {
                                        // if input file was prepared
					// under DOS/Windows
                                        In.get(cIn);
                                        if (In.eof()) {
                                                return -1;
                                        }
                                }
                        }
                }
                goto skip_again;

        case '/':
                In.get(cIn);
                if (In.eof()) {
                        return -1;
		}

                if (cIn == '*') {
                        // ok, it's the beginning of a C-style comment, '/*'; eat until '*/'
			for (;;) {
				In.get(cIn);
				if (In.eof()) {
					return -1;
				}
                                if (cIn == '*') {
end_of_comment:;
                                        // there may be more than just one '*'
                                        // before the end of a comment...
                                        while (cIn == '*') {
                                                In.get(cIn);
                                                if (In.eof()) {
                                                        return -1;
                                                }
                                                if (cIn == '/') {
							// C-style comment is over
                                                        goto skip_again;
                                                }
                                        }

                                } else if (cIn == '/') {
                                        In.get(cIn);
                                        if (In.eof()) {
                                                return -1;
                                        }
                                        if (cIn == '*') {
                                                silent_cerr("warning: '/*' inside a comment "
                                                        "at line " << HP.GetLineData()
                                                        << std::endl);
                                                goto end_of_comment;
                                        }
                                }
                        }

                } else {
                        In.putback(cIn);
                        return 0;
                }
        }

        return 0;
}

LowParser::LowParser(HighParser& hp)
: HP(hp)
{
        NO_OP;
}

LowParser::~LowParser(void)
{
        NO_OP;
}

void
LowParser::PackWords(InputStream& In)
{
       	char cIn;

        sCurrWord.clear();

        /* note: no remarks allowed inside words */
        for (;;) {
		In.get(cIn);
		if (In.eof()) {
			return;
		}
                switch (cIn) {
                case COLON:
                case COMMA:
                case SEMICOLON:
                        goto end_of_word;

                default:
                        if (!isspace(cIn)) {
                                sCurrWord += cIn;
                        }
                }
        }

        throw EndOfFile(MBDYN_EXCEPT_ARGS);

end_of_word:;

        In.putback(cIn);
}

LowParser::Token
LowParser::GetToken(InputStream& In)
{
        /* toglie gli spazi iniziali e tutti i commenti */
        char cIn;
        if (skip_remarks(HP, In, cIn)) {
                return CurrToken = LowParser::ENDOFFILE;
        }

        if (isalpha(cIn) || cIn == '_') {
                PackWords(In.putback(cIn));
                return CurrToken = LowParser::WORD;
        }

        switch (cIn) {
        case ',':
                return CurrToken = LowParser::COMMA;

        case ':':
                return CurrToken = LowParser::COLON;

        case ';':
                return CurrToken = LowParser::SEMICOLON;

        case '.':
        case '-':
        case '+':
is_digit:;
                In.putback(cIn) >> dCurrNumber;
                return CurrToken = LowParser::NUMBER;

        default:
                if (isdigit(cIn)) {
                        goto is_digit;
                }
                In.putback(cIn);
                return CurrToken = LowParser::UNKNOWN;
        }
}


doublereal
LowParser::dGetReal(void) const
{
        return dCurrNumber;
}


integer
LowParser::iGetInt(void) const
{
        return integer(dCurrNumber);
}


const char*
LowParser::sGetWord(void)
{
        return sCurrWord.c_str();
}

/* LowParser - end */


/* KeyTable - begin */

KeyTable::KeyTable(HighParser& hp, const char* const sTable[])
: sKeyWords(0), oldKey(0), HP(hp)
{
        sKeyWords = (char* const*)sTable;
        oldKey = HP.PutKeyTable(*this);
}


KeyTable::~KeyTable(void)
{
        if (oldKey) {
                (void)HP.PutKeyTable(*oldKey);
        }
}


int
KeyTable::Find(const char* sToFind) const
{
        for (int iCnt = 0; sKeyWords[iCnt]; iCnt++) {
                if (strcasecmp(sKeyWords[iCnt], sToFind) == 0) {
                        return iCnt;
                }
        }

        return -1;
}

const char* KeyTable::pGetDescription(integer iIndex) const
{
#ifdef DEBUG
     for (integer i = 0; i < iIndex; ++i) {
          ASSERT(sKeyWords[i]);
     }
#endif
     return sKeyWords[iIndex];
}

/* KeyTable - end */


/* DescRead - begin */

/* bag that contains functions to parse descriptions */

typedef std::map<std::string, DescRead *, ltstrcase> DescFuncMapType;
static DescFuncMapType DescFuncMap;

struct DescWordSetType : public HighParser::WordSet {
        bool IsWord(const std::string& s) const {
                return DescFuncMap.find(s) != DescFuncMap.end();
        };
};

static DescWordSetType DescWordSet;

bool
SetDescData(const std::string& name, DescRead *rf)
{
        pedantic_cout("registering description \"" << name << "\"" << std::endl);
        return DescFuncMap.insert(DescFuncMapType::value_type(name, rf)).second;
}

/* Reads descriptions */

bool
ReadDescription(HighParser& HP, const std::string& desc)
{
        DEBUGCOUTFNAME("ReadDescription()");

        bool bRC(false);
        DescFuncMapType::iterator func = DescFuncMap.find(desc);
        if (func != DescFuncMap.end()) {
                HP.GotDescription();
                if (!HP.IsArg() && !HP.IsDescription()) {
                        silent_cerr("Parser error in ReadDescription(),"
                                " colon or semicolon expected after description at line "
                                << HP.GetLineData() << std::endl);
                        throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
                }

                bRC = func->second->Read(HP);

                if (HP.IsArg()) {
                        silent_cerr("semicolon expected at line " << HP.GetLineData() << std::endl);
                        throw HighParser::ErrSemicolonExpected(MBDYN_EXCEPT_ARGS);
                }
        }

        return bRC;
}

DescRead::~DescRead(void)
{
        NO_OP;
}

struct RemarkDR : public DescRead {
public:
        bool Read(HighParser& HP);
};

bool
RemarkDR::Read(HighParser& HP)
{
        silent_cout("line " << HP.GetLineData());

        char prefix = ':';
        while (HP.IsArg()) {
                TypedValue v;
                v = HP.GetValue(v);
                silent_cout(prefix << ' ' << v);

                if (prefix == ':') {
                        prefix = ',';
                }
        }

        silent_cout(std::endl);

        return true;
}

struct PrintSymbolTableDR : public DescRead {
public:
        bool Read(HighParser& HP);
};

bool
PrintSymbolTableDR::Read(HighParser& HP)
{
        if (!HP.IsArg()) {
                silent_cout( "math parser symbol table at line "
                        << HP.GetLineData() << ":" << std::endl
                        << HP.GetMathParser().GetSymbolTable() << std::endl);
                return true;
        }

        if (HP.IsKeyWord("all")) {
                const MathParser::NameSpaceMap& ns = HP.GetMathParser().GetNameSpaceMap();
                for (MathParser::NameSpaceMap::const_iterator i = ns.begin(); i != ns.end(); ++i) {
                        const std::string& sName = i->second->sGetName();
                        const Table *pT = i->second->GetTable();
                        if (pT != 0) {
                                silent_cout( "namespace \"" << sName << "\" symbol table at line "
                                        << HP.GetLineData() << ":" << std::endl
                                        << *pT << std::endl);
                        }
                }

                return true;
        }

        while (HP.IsArg()) {
                const std::string sName = HP.GetString();
                MathParser::NameSpace *pN = HP.GetMathParser().GetNameSpace(sName);
                if (pN == 0) {
                        silent_cerr("PrintSymbolTableDR::Read(): warning, unable to find namespace \"" << sName << "\" at line "
                                << HP.GetLineData() << std::endl);

                } else {
                        Table *pT = pN->GetTable();
                        if (pT == 0) {
                                silent_cerr("PrintSymbolTableDR::Read(): warning, namespace \"" << sName << "\" "
                                        "has no symbol table at line " << HP.GetLineData() << std::endl);

                        } else {
                                silent_cout( "namespace \"" << sName << "\" symbol table at line "
                                        << HP.GetLineData() << ":" << std::endl
                                        << *pT << std::endl);
                        }
                }
        }

        return true;
}

struct SetDR : public DescRead {
public:
        bool Read(HighParser& HP);
};

bool
SetDR::Read(HighParser& HP)
{
        if (!HP.IsArg()) {
                silent_cerr("Parser error in SetDR::Read(), "
                        "arg expected at line "
                        << HP.GetLineData() << std::endl);
                throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
        }

        TypedValue v;
        HP.GetValue(v);

        return true;
}

struct SetEnvDR : public DescRead {
public:
        bool Read(HighParser& HP);
};

bool
SetEnvDR::Read(HighParser& HP)
{
#ifdef HAVE_SETENV
        if (!HP.IsArg()) {
                silent_cerr("Parser error in SetEnvDR::Read(), "
                        "arg(s) expected at line "
                        << HP.GetLineData() << std::endl);
                throw HighParser::ErrColonExpected(MBDYN_EXCEPT_ARGS);
        }

        int overwrite = 0;
        if (HP.IsKeyWord("overwrite")) {
                bool b = HP.GetYesNoOrBool();
                overwrite = b ? 1 : 0;
        }

        const std::string ava = HP.GetStringWithDelims();
        if (ava.empty()) {
                silent_cerr("unable to get AVA for \"setenv\" at line "
                                << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const std::string::size_type sep = ava.find('=');
        if (sep == std::string::npos) {
#ifdef HAVE_UNSETENV
                unsetenv(ava.c_str());
#elif defined(HAVE_PUTENV)
		static std::list<std::string> environmentStorage;
		environmentStorage.push_back(ava);
		if (putenv(environmentStorage.back().data())) {
			environmentStorage.pop_back();
                        silent_cerr("unable to unset the environment variable "
                                        "\"" << ava << "\" at line "
                                        << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
#endif	/* !HAVE_UNSETENV && !HAVE_PUTENV */

        } else {
                if (sep == 0) {
                        silent_cerr("illegal AVA \"" << ava
                                        << "\" at line "
                                        << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

		const std::string name = ava.substr(0, sep);
		const std::string value = ava.substr(sep + 1);
		bool bPresent(getenv(name.c_str()) != NULL);
		int rc = setenv(name.c_str(), value.c_str(), overwrite);
                if (rc) {
                        silent_cerr("unable to set the environment variable \""
					<< name << "\" to \"" << value
                                        << "\" at line " << HP.GetLineData()
                                        << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (bPresent && overwrite == 0) {
			silent_cout("Environment variable \"" << name
				<< "\" _not_ overwritten with \"" << value
				<< "\" (current value is \"" << getenv(name.c_str())
                                << "\") at line " << HP.GetLineData()
                                << std::endl);

                } else if (!bPresent) {
			silent_cout("Environment variable \"" << name
				<< "\" set to \"" << value
                                << "\" at line " << HP.GetLineData()
                                << std::endl);

                } else {
			silent_cout("Environment variable \"" << name
				<< "\" overwritten to \"" << value
                                << "\" at line " << HP.GetLineData()
                                << std::endl);
                }
        }
#else // ! HAVE_SETENV
        silent_cerr("SetEnvDR::Read(): warning, setenv() not available; ignored at line "
                << HP.GetLineData() << std::endl);
#endif // !HAVE_SETENV
        return true;
}

struct ExitDR : public DescRead {
public:
        bool Read(HighParser& HP);
};

bool
ExitDR::Read(HighParser& HP)
{
        if (!HP.IsDescription()) {
                silent_cerr("Parser error in ExitDR::Read(),"
                        " semicolon expected at line "
                        << HP.GetLineData() << std::endl);
                throw HighParser::ErrSemicolonExpected(MBDYN_EXCEPT_ARGS);
        }

        /* exits with no error */
        throw NoErr(MBDYN_EXCEPT_ARGS);
}

static unsigned desc_done;

static void
InitDescData(void)
{
        if (::desc_done++ > 0) {
                return;
        }

        SetDescData("remark", new RemarkDR);
        SetDescData("print" "symbol" "table", new PrintSymbolTableDR);
        SetDescData("set", new SetDR);
        SetDescData("setenv", new SetEnvDR);
        SetDescData("exit", new ExitDR);

        /* NOTE: add here initialization of new built-in descriptions;
         * alternative ways to register new custom descriptions are:
         * - call SetDescData() from anywhere in the code
         * - write a module that calls SetDescData() from inside a function
         *   called module_init(), and run-time load it using "module load"
         *   in the input file.
         */
}

static void
DestroyDescData(void)
{
        if (::desc_done == 0) {
                silent_cerr("DestroyDescData() called once too many" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (--::desc_done > 0) {
                return;
        }

        /* free stuff */
        for (DescFuncMapType::iterator i = DescFuncMap.begin(); i != DescFuncMap.end(); ++i) {
                delete i->second;
        }
        DescFuncMap.clear();
}

/* DescRead - end */


/* HighParser - begin */

static std::stack<const HighParser *> pHP;

static const HighParser::ErrOut unknownErr = { "(unknown)", "(unknown)", 0 };

HighParser::ErrOut
mbdyn_get_line_data(void)
{
        if (!pHP.empty()) {
                return pHP.top()->GetLineData();
        }

        return unknownErr;
}

std::ostream&
mbdyn_print_line_data(std::ostream& out)
{
        if (!pHP.empty()) {
                out << pHP.top()->GetLineData();
        }

        return out;
}

HighParser::HighParser(MathParser& MP, InputStream& streamIn)
: ESCAPE_CHAR('\\'),
LowP(*this),
pIn(&streamIn),
pf(NULL),
MathP(MP),
KeyT(0)
{
        DEBUGCOUTFNAME("HighParser::HighParser");
        CurrToken = HighParser::DESCRIPTION;

        InitDescData();

        pHP.push(this);
}


HighParser::~HighParser(void)
{
        DEBUGCOUTFNAME("HighParser::~HighParser");
        Close();
        ASSERT(pHP.top() == this);
        pHP.pop();

        DestroyDescData();
}


void
HighParser::Close(void)
{
        NO_OP;
}


const KeyTable*
HighParser::PutKeyTable(const KeyTable& KT)
{
        const KeyTable* oldKey = KeyT;

        KeyT = &KT;

        return oldKey;
}

MathParser&
HighParser::GetMathParser(void)
{
        return MathP;
}

int
HighParser::GetLineNumber(void) const
{
        return const_cast<InputStream *>(pIn)->GetLineNumber();
}


HighParser::ErrOut
HighParser::GetLineData(void) const
{
        ErrOut LineData;
        LineData.iLineNumber = GetLineNumber();
        LineData.sFileName.clear();
        LineData.sPathName.clear();
        return LineData;
}


bool
HighParser::IsDescription(void) const
{
        return (CurrToken == HighParser::DESCRIPTION);
}

HighParser::Token
HighParser::GotDescription(void)
{
        return FirstToken();
}

int
HighParser::iGetDescription_int(const char* const s)
{
        int i = -1;

        if (KeyT) {
                i = KeyT->Find(s);
        }

        if (FirstToken() == HighParser::UNKNOWN) {
                silent_cerr("Parser error in HighParser::iGetDescription_int(), "
                        "semicolon expected at line "
                        << GetLineData() << std::endl);
                throw HighParser::ErrSemicolonExpected(MBDYN_EXCEPT_ARGS);
        }

        return i;
}


void
HighParser::Eof(void)
{
         throw EndOfFile(MBDYN_EXCEPT_ARGS);
}

int
HighParser::GetDescription(void)
{
        /* Checks if current token is a description */
        if (!IsDescription()) {
                silent_cerr("Parser error in HighParser::GetDescription, "
                        "invalid call to GetDescription at line "
                        << GetLineData() << std::endl);
                throw HighParser::ErrInvalidCallToGetDescription(MBDYN_EXCEPT_ARGS);
        }

restart_parsing:;

        CurrLowToken = LowP.GetToken(*pIn);
        if (CurrLowToken != LowParser::WORD) {
                if (pIn->eof()) {
                        Eof();
                        goto restart_parsing;
                }

                silent_cerr("Parser error in HighParser::GetDescription, "
                        << "keyword expected at line "
                        << GetLineData() << std::endl);
                throw HighParser::ErrKeyWordExpected(MBDYN_EXCEPT_ARGS);
        }

        /* Description corrente */
        const char* s = LowP.sGetWord();

        if (ReadDescription(*this, s)) {
                goto restart_parsing;
        }

        return iGetDescription_int(s);
}


HighParser::Token
HighParser::FirstToken(void)
{
        CurrLowToken = LowP.GetToken(*pIn);

        switch (CurrLowToken) {
        case LowParser::COLON:
                CurrToken = HighParser::ARG;
                break;

        case LowParser::SEMICOLON:
                CurrToken = HighParser::DESCRIPTION;
                break;

        default:
                CurrToken = HighParser::UNKNOWN;
                break;
        }

        return CurrToken;
}

void
HighParser::ExpectDescription(void)
{
        /* forces the next expected token to be a "description"
         * e.g. a keyword followed by a colon (deprecated) */
        CurrToken = HighParser::DESCRIPTION;
}


void
HighParser::ExpectArg(void)
{
        /* forces the next expected token to be an argument
         * e.g. a keyword followed by a separator (deprecated) */
        CurrToken = HighParser::ARG;
}


bool
HighParser::IsArg(void)
{
        return (CurrToken == ARG);
}

void
HighParser::PutBackSemicolon(void)
{
        if (CurrLowToken == LowParser::SEMICOLON) {
                pIn->putback(';');
        }
}


void
HighParser::NextToken(const char* sFuncName)
{
        CurrLowToken = LowP.GetToken(*pIn);
        switch (CurrLowToken) {
        case LowParser::COMMA:
                CurrToken = HighParser::ARG;
                break;

        case LowParser::SEMICOLON:
                CurrToken = HighParser::DESCRIPTION;
                break;

        default:
                silent_cerr("Parser error in "
                        << sFuncName << ", missing separator at line "
                        << GetLineData() << std::endl);
                throw HighParser::ErrMissingSeparator(MBDYN_EXCEPT_ARGS);
        }
}

int
HighParser::ParseWord(unsigned flags)
{
        sStringBuf.clear();
        sStringBufWithSpaces.clear();

        char cIn;
        if (skip_remarks(*this, *pIn, cIn)) {
                return CurrToken = HighParser::ENDOFFILE;
        }

        if (!isalpha(cIn) && cIn != '_') {
                pIn->putback(cIn);
                return -1;
        }

        sStringBufWithSpaces += cIn;

        if (flags & LOWER) {
                sStringBuf += tolower(cIn);

        } else if (flags & UPPER) {
                sStringBuf += toupper(cIn);

        } else {
                sStringBuf += cIn;
        }

	for (;;) {
		pIn->get(cIn);
		if (pIn->eof()) {
			return -1;
		}

		if (!(isalnum(cIn) || cIn == '_' || isspace(cIn))) {
			break;
		}

                sStringBufWithSpaces += cIn;

                if (isspace(cIn)) {
                        continue;
                }

		char c = cIn;
                if (flags & LOWER) {
                        c = tolower(c);

                } else if (flags & UPPER) {
                        c = toupper(c);
		}

                sStringBuf += c;
        }

        pIn->putback(cIn);

        return 0;
}

void
HighParser::PutbackWord(void)
{
        for (std::string::reverse_iterator i = sStringBufWithSpaces.rbegin();
                i != sStringBufWithSpaces.rend(); ++i)
        {
                pIn->putback(*i);
        }
}

bool
HighParser::IsKeyWord(const char* sKeyWord)
{
        const char sFuncName[] = "HighParser::IsKeyWord()";

        if (CurrToken != HighParser::ARG) {
                return false;
        }

        switch (ParseWord()) {
        case 0:
                break;

        case HighParser::ENDOFFILE:
                return true;

        default:
                return false;
        }

        if (!strcasecmp(sStringBuf.c_str(), sKeyWord)) {
                NextToken(sFuncName);
                return true;
        }

        PutbackWord();

        return false;
}

int
HighParser::IsKeyWord(void)
{
        const char sFuncName[] = "HighParser::IsKeyWord()";

        if (CurrToken != HighParser::ARG) {
                return -1;
        }

        switch (ParseWord()) {
        case 0:
                break;

        case HighParser::ENDOFFILE:
                return HighParser::ENDOFFILE;

        default:
                return -1;
        }

        int iKW = -1;

        if (KeyT) {
                iKW = KeyT->Find(sStringBuf.c_str());
        }

        if (iKW >= 0) {
                NextToken(sFuncName);
                return iKW;
        }

        PutbackWord();

        return -1;
}

/* 1 se l'argomento successivo e' una parola in un WordSet */
const char *
HighParser::IsWord(const HighParser::WordSet& ws)
{
        const char sFuncName[] = "HighParser::IsWord()";

        if (CurrToken != HighParser::ARG) {
                return 0;
        }

        switch (ParseWord()) {
        case 0:
                break;

        default:
                return 0;
        }

        if (ws.IsWord(sStringBuf)) {
                NextToken(sFuncName);
                return sStringBuf.c_str();
        }

        PutbackWord();

        return 0;
}

TypedValue
HighParser::GetValue(const TypedValue& vDefVal)
{
        return GetValue<range_any<TypedValue> >(vDefVal, range_any<TypedValue>());
}

bool
HighParser::GetBool(bool bDefVal)
{
        TypedValue v(bDefVal);
        v = GetValue(v);
        return v.GetBool();
}

/*
 * Read a boolean as "yes" or "no" and put the result in bRet
 * return true in case of success, false otherwise (bRet undefined)
 */
bool
HighParser::GetYesNo(bool& bRet)
{
        if (IsKeyWord("yes")) {
                bRet = true;

        } else if (IsKeyWord("no")) {
                bRet = false;

        } else {
                return false;
        }

        return true;
}

bool
HighParser::GetYesNoOrBool(bool bDefval)
{
        bool bRet;

        if (!GetYesNo(bRet)) {
                bRet = GetBool(bDefval);
        }

        return bRet;
}

integer
HighParser::GetInt(integer iDefVal)
{
        return GetInt<range_any<integer> >(iDefVal, range_any<integer>());
}

doublereal
HighParser::GetReal(const doublereal& dDefVal)
{
        return GetReal<range_any<doublereal> >(dDefVal, range_any<doublereal>());
}

mbsleep_t
HighParser::GetTimeout(const mbsleep_t& DefVal)
{
        doublereal d;
        mbsleep_sleep2real(&DefVal, &d);
        TypedValue v(d);
        v = GetValue(v);
        mbsleep_t newval;
        mbsleep_real2sleep(v.GetReal(), &newval);
        return newval;
}

std::string
HighParser::GetString(const std::string& sDefVal)
{
        TypedValue v(sDefVal);
        v = GetValue(v);
        return v.GetString();
}


int
HighParser::GetWord(void)
{
        const char sFuncName[] = "HighParser::GetWord()";

        if (CurrToken != HighParser::ARG) {
                silent_cerr("Parser error in "
                        << sFuncName << ", keyword arg expected at line "
                        << GetLineData() << std::endl);
                throw HighParser::ErrKeyWordExpected(MBDYN_EXCEPT_ARGS);
        }

        CurrLowToken = LowP.GetToken(*pIn);
        if (CurrLowToken != LowParser::WORD) {
                silent_cerr("Parser error in "
                        << sFuncName << ", keyword expected at line "
                        << GetLineData() << std::endl);
                throw HighParser::ErrKeyWordExpected(MBDYN_EXCEPT_ARGS);
        }

        int i = -1;
        if (KeyT) {
                i = KeyT->Find(LowP.sGetWord());
        }

        NextToken(sFuncName);

        return i;
}

std::string
HighParser::GetString(unsigned flags)
{
        const char sFuncName[] = "HighParser::GetString()";

        pedantic_cout("use of deprecated method \"GetString\" at line"
                << GetLineData() << std::endl);

        if (CurrToken != HighParser::ARG) {
                silent_cerr("Parser error in "
                        << sFuncName << ", string arg expected at line "
                        << GetLineData() << std::endl);
                throw HighParser::ErrStringExpected(MBDYN_EXCEPT_ARGS);
        }

        std::string result;

        char cIn = '\0';
	for (;;) {
		pIn->get(cIn);
		if (pIn->eof()) {
                	CurrToken = HighParser::ENDOFFILE;
			return result;
		}

		if (!isspace(cIn)) {
			break;
		}
	}

        pIn->putback(cIn);
	for (;;) {
		pIn->get(cIn);
		if (pIn->eof()) {
                        CurrToken = HighParser::ENDOFFILE;
                        return result;
		}

		if (!(cIn != ',' && cIn != ';')) {
			break;
		}

                if (!(flags & HighParser::EATSPACES) || !isspace(cIn)) {
			char c = cIn;
                        if (flags & HighParser::LOWER) {
                                c = tolower(c);

                        } else if (flags & HighParser::UPPER) {
                                c = toupper(c);
                        }

                        result += c;
                }
        }

        pIn->putback(cIn);

        NextToken(sFuncName);

        return result;
}

void
HighParser::SetDelims(enum Delims Del, char &cLdelim, char &cRdelim) const
{
        cLdelim = '\0';
        cRdelim = '\0';

        switch (Del) {
        case PLAINBRACKETS:
                cLdelim = '(';
                cRdelim = ')';
                break;

        case SQUAREBRACKETS:
                cLdelim = '[';
                cRdelim = ']';
                break;

        case CURLYBRACKETS:
                cLdelim = '{';
                cRdelim = '}';
                break;

        case SINGLEQUOTE:
                cLdelim = '`';
                cRdelim = '\'';
                break;

        default:
        case UNKNOWNDELIM:
        case DEFAULTDELIM:
        case DOUBLEQUOTE:
                cLdelim = '"';
                cRdelim = '"';
                break;
        }
}

bool
HighParser::IsStringWithDelims(enum Delims Del)
{
        char cLdelim, cRdelim;
        SetDelims(Del, cLdelim, cRdelim);

        char cIn;
        if (skip_remarks(*this, *pIn, cIn)) {
                return false;
        }

        /* put back the first non-remark char */
        pIn->putback(cIn);

        /* if the left delimiter is found, true */
        return (cIn == cLdelim);
}

std::string
HighParser::GetStringWithDelims(enum Delims Del, bool escape)
{
        const char sFuncName[] = "HighParser::GetStringWithDelims()";

        if (CurrToken != HighParser::ARG) {
                silent_cerr("Parser error in "
                        << sFuncName << ", string arg expected at line "
                        << GetLineData() << std::endl);
                throw HighParser::ErrStringExpected(MBDYN_EXCEPT_ARGS);
        }

        std::string result;

        char cLdelim, cRdelim;
        SetDelims(Del, cLdelim, cRdelim);

        char cIn;
        if (skip_remarks(*this, *pIn, cIn)) {
                return result;
        }

        // Se trova il delimitatore sinistro, legge la stringa
        if (cIn == cLdelim) {
                for (;;) {
			pIn->get(cIn);
                        if (pIn->eof()) {
                                // FIXME: this should be an error ...
				silent_cerr("End-of-file encountered in " << sFuncName << " while looking for right string delimiter '" << cRdelim << "' after " << unsigned(result.size()) << " characters at line " << GetLineData() << std::endl);
                		throw EndOfFile(MBDYN_EXCEPT_ARGS);
			}

			if (cIn == cRdelim) {
				break;
			}

                        if (cIn == ESCAPE_CHAR) {
                                pIn->get(cIn);
                		if (pIn->eof()) {
                		        // FIXME: this should be an error ...
					silent_cerr("End-of-file encountered in " << sFuncName << " while escaping a char after " << unsigned(result.size()) << " characters at line " << GetLineData() << std::endl);
        				throw EndOfFile(MBDYN_EXCEPT_ARGS);
				}

                                if (cIn == '\n') {

                                        /*
                                         * eat the newline as well, so that

                                                "first line\
                                                second line"

                                         * actually results in "first linesecond line"
                                         */

                                        pIn->get(cIn);
        				if (pIn->eof()) {
        		       			// FIXME: this should be an error ...
						silent_cerr("End-of-file encountered in " << sFuncName << " while escaping a newline ('\\n') after " << unsigned(result.size()) << " characters at line " << GetLineData() << std::endl);
        					throw EndOfFile(MBDYN_EXCEPT_ARGS);
					}

                                } else if (cIn == '\r') {
                                        pIn->get(cIn);
        				if (pIn->eof()) {
        		       			// FIXME: this should be an error ...
						silent_cerr("End-of-file encountered in " << sFuncName << " while escaping a 'carriage return' ('\\r') after " << unsigned(result.size()) << " characters at line " << GetLineData() << std::endl);
        					throw EndOfFile(MBDYN_EXCEPT_ARGS);
					}

                                        if (cIn != '\n') {
                                                pIn->putback(cIn);
                                                goto escaped_generic;
                                        }

                                        pIn->get(cIn);
        				if (pIn->eof()) {
        		       			// FIXME: this should be an error ...
						silent_cerr("End-of-file encountered in " << sFuncName << " after escaping a 'carriage return' ('\\r') after " << unsigned(result.size()) << " characters at line " << GetLineData() << std::endl);
        					throw EndOfFile(MBDYN_EXCEPT_ARGS);
					}

                                } else if ((cIn == ESCAPE_CHAR) || (cIn == cRdelim)) {
                                        if (!escape) {
                                                result += ESCAPE_CHAR;
                                        }

                                } else {
escaped_generic:;
                                        if (escape) {
                                                int i, c = 0;
                                                char hex[3];

                                                /*
                                                 * allow non-printable chars in the form "\<hexpair>",
                                                 * so that "\78" is equivalent to "x";
                                                 * "\<non-hexpair>" is treated as an error.
                                                 */

                                                hex[0] = cIn;
						pIn->get(cIn);
        					if (pIn->eof()) {
        		       				// FIXME: this should be an error ...
							silent_cerr("End-of-file encountered in " << sFuncName << " while escaping a hexpair char ('\\dd', where 'dd' are two hex digits) after " << unsigned(result.size()) << " characters at line " << GetLineData() << std::endl);
        						throw EndOfFile(MBDYN_EXCEPT_ARGS);
						}
                                                hex[1] = cIn;
                                                hex[2] = '\0';

                                                for (i = 0; i < 2; i++) {
                                                        int shift = 4*(1 - i), h = 0;

                                                        /* NOTE: this conversion relies
                                                         * on 0-9, a-f, A-F being consecutive,
                                                         * which is true for ASCII, but might
                                                         * not be for other encodings;
                                                         * bah, not critical right now */
                                                        if (hex[i] >= '0' && hex[i]  <= '9') {
                                                                h = hex[i] - '0';
                                                        } else if (hex[i] >= 'a' && hex[i] <= 'f') {
                                                                h = hex[i] - 'a';
                                                        } else if (hex[i] >= 'A' && hex[i] <= 'F') {
                                                                h = hex[i] - 'A';
                                                        } else {
                                                                silent_cerr("invalid escape sequence "
                                                                        "\"\\" << hex << "\" "
                                                                        "at line " << GetLineData()
                                                                        << std::endl);
                                                                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                                        }

                                                        c += (h << shift);
                                                }
                                                cIn = c;

                                        } else {
                                                result += ESCAPE_CHAR;
                                        }
                                }
                        }
                        result += cIn;
                }

                /* Se trova una virgola o un punto e virgola, le rimette nello stream
                 * e passa oltre, restituendo una stringa vuota. Il chiamante deve
                 * occuparsi della gestione del valore di default */
        } else if (cIn == ',' || cIn == ';') {
                pIn->putback(cIn);

                /* Altrimenti c'e' qualcosa senza delimitatore. Adesso da' errore,
                 * forse e' piu' corretto fargli ritornare lo stream intatto */
        } else {
                silent_cerr("Parser error in "
                        << sFuncName << std::endl
                        << "first non-blank char at line "
                        << GetLineData() << " isn't a valid left-delimiter"
                        << std::endl);
                throw HighParser::ErrIllegalDelimiter(MBDYN_EXCEPT_ARGS);
        }

        NextToken(sFuncName);
        return result;
}

/* Returns the current input stream */
InputStream&
HighParser::GetInputStream(void) const
{
        ASSERT(pIn != 0);
        return *pIn;
}

/* Sets a new input stream (use with care!) */
const void
HighParser::PutInputStream(InputStream& streamIn)
{
        pIn = &streamIn;
}

/* HighParser - end */

std::ostream&
operator << (std::ostream& out, const HighParser::ErrOut& err)
{
        out << err.iLineNumber;

        if (!err.sFileName.empty()) {
                out << ", file <";
                if (!err.sPathName.empty()) {
                        out << err.sPathName << DIR_SEP;
                }
                out << err.sFileName << '>';
        }

        return out;
}

/* HighParser - end */
