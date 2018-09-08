#ifndef __TE_H
#define __TE_H

#include <fstream>
#include <iostream>     // std::cout
#include <algorithm>
#include <regex>
#include <map>
#include <utility>      // std::pair

#include <string>       // std::string
#include <streambuf>
#include <sstream>      // std::ostringstream

#ifndef NDEBUG
#define  DEBUG
#endif


using TRenderItem = std::map<std::string, std::string>;
using TRenderData = std::multimap<std::string, TRenderItem>;

using namespace std::string_literals; // @suppress("Using directive in header file")

class Cte
    {
	using TMapOfBlocks = std::map<std::string, std::string>;
	using TLoopParams  = std::pair<std::string, std::string>;

	std::string const m_sFilepath;
	std::string m_sPage;

    public:

	Cte(TRenderData & mData, std::string const & crsFilename, std::string const & crsFilepath = "")
	    : m_sFilepath(crsFilepath)
	    {
	    std::string sTemplate = ReadTemplate(m_sFilepath + crsFilename);

	    sTemplate = GetIfs(sTemplate, mData);
            sTemplate = FillLoops(sTemplate, mData);
	    sTemplate = FillVariables(sTemplate, mData);

	    std::smatch      sm{};
	    std::regex const re(R"e(\{\%\s+extends\s+\"(.*)\"\s+\%\})e");
	    std::string      sExtend{};
	    std::regex_search(sTemplate, sm, re);
	    if ( sm.size() > 1 )
		{
		sExtend = sm[1];

		TMapOfBlocks const mBlocks = GetBlocks(sTemplate);

		sTemplate = ReadTemplate(crsFilepath + sExtend);
	        sTemplate = GetIfs(sTemplate, mData);
                sTemplate = FillLoops(sTemplate, mData);
		sTemplate = FillVariables(sTemplate, mData);
		sTemplate = FillBlocks(sTemplate, mBlocks);
		}
	    m_sPage = std::move(sTemplate);
	    }

	size_t length() const { return m_sPage.length(); }

	/**
	 * @brief The output stream operator for Cte outputs the resulting page
	 *
	 * @param ros The output stream
	 * @param crote The template engine
	 */
	friend std::ostream & operator << (std::ostream & ros, Cte const & crTE)
	    {
	    return ros << crTE.m_sPage;
	    }

	/**
	 * @brief Converts " value in values %}" to a pair<"value", "values">
	 * 	  where values is the key of the values, who go into {{ value }}
	 *
	 * @param crsFragment The statement fragment to be sorted
	 */
        TLoopParams SplitLoopParam(std::string const & crsFragment) const
            {

	    std::smatch      sm{};
	    std::regex const re(R"(([^\s*]*)\s+in\s+([^\s*]*)\s+\%\})");
	    std::regex_search(crsFragment, sm, re);
	    if ( sm.size() > 2 )
		{
		return std::make_pair(sm[1], sm[2]);
		}
            return std::make_pair("", "");
            }

        /**
         * @brief Replaces "{% for message in messages %}<li>{{ message }}</li>{% endfor %}"
         *        with "<li>`mVListsat("message")`</li>" for each message in messages
         */
        std::string FillLoops(std::string const & crsPage, TRenderData const & mm) const
	    {
//	    { {"property", { {"id", "1"}, {"", "nm9087684"} } },
//	      {"property", { {"id", "2"}, {"", "actor"}     } } }
//
//          "{% for property i dummy %}{{ property.id }}, {% endfor %}" => "1, 2, "

	    std::ostringstream oss{};

	    std::regex const re(R"(\{\%\s+for\s+|\{\%\s+endfor\s+\%\})");
	    size_t n{0};
	    for (auto it = std::sregex_token_iterator(crsPage.begin(), crsPage.end(), re, -1); it != std::sregex_token_iterator(); ++it)
		{
		if ( ++n & 1 )
		    {
		    oss << *it;
		    }
		else
		    {
		    std::string const si    = *it;
		    size_t      const p     = si.find('}')+1;
                    auto        const param = SplitLoopParam( si.substr(0, p) );
                    auto        const er    = mm.equal_range(param.first);
#ifdef DEBUG
                   // how to detect an error
#endif
                    for (auto itm=er.first; itm!=er.second; ++itm)
                	{
                        oss << FillVariables(si.substr(p), param.first, itm->second);
                	}
		    }
		}
	    return oss.str();
	    }

        /**
         * @brief Replaces '{{ crsName }}' with 'value'
         *
         * Deals with repeatable sequences from between enclosing tags like
         * for-loops
         *
         * @param crsPart A part of the document
         * @param crsName The variables base name
         * @param mVariables The data
         */
        std::string FillVariables(std::string const & crsPart,
        			  std::string const & crsName,
				  TRenderItem     const & mVariables) const
            {
	    std::ostringstream oss{};

	    std::regex const re(R"(\{{2}\s+|\s+\}{2})");
	    size_t n{0};
	    for (auto it = std::sregex_token_iterator(crsPart.begin(), crsPart.end(), re, -1); it != std::sregex_token_iterator(); ++it)
		{
		if ( ++n & 1 )
		    {
		    oss << *it;
		    }
		else
		    {
//                  {"my", { {"", "demo"},    {"pk",   "demo-private"} } },
//		 =>     {{ my }} = demo & {{ my.pk }} = demo-private

		    std::string sk1 = *it;
		    std::string sk2 = "";
		    size_t const  p = sk1.find('.');
		    if ( p != std::string::npos )
			{
			sk2 = sk1.substr(p+1);
			sk1 = sk1.substr(0,p);
			}

		    if ( crsName != sk1 )
			{
			oss << "{{ !!!" << *it << "!!! }}";
			continue;
			}

		    auto const a = mVariables.find(sk2);
		    if ( a != mVariables.end() )
			{
			oss << a->second;
			}
#ifdef DEBUG
                    else { std::cerr << "\nte-debug, variable [" << sk2 << "] not found\n"; }
#endif
		    }
		}
	    return oss.str();
	    }

	/**
	 * @brief Replaces "{{ value }}" with the value.
	 *
	 * @param crsPage The template
	 * @param mData The data to fill in
	 */
        std::string FillVariables(std::string const & crsPage, TRenderData& mData) const
	    {
	    std::ostringstream oss{};

	    std::regex const re(R"(\{{2}\s+|\s+\}{2})");
	    size_t n{0};
	    for (auto it = std::sregex_token_iterator(crsPage.begin(), crsPage.end(), re, -1); it != std::sregex_token_iterator(); ++it)
		{
		if ( ++n & 1 )
		    {
		    oss << *it;
		    }
		else
		    {
//                  {"my", { {"", "demo"},    {"pk",   "demo-private"} } },
//		 =>     {{ my }} = demo & {{ my.pk }} = demo-private
//		    oss << mVariables.find(sk1)->second.find(sk2)->second;

		    std::string sk1 = *it;
		    std::string sk2 = "";
		    size_t const  p = sk1.find('.');
		    if ( p != std::string::npos )
			{
			sk2 = sk1.substr(p+1);
			sk1 = sk1.substr(0,p);
			}

		    auto const a = mData.find(sk1);
		    if ( a != mData.end() )
			{
			auto const b = a->second.find(sk2);
			if ( b != a->second.end() )
			    {
			    oss << b->second;
			    }
#ifdef DEBUG
                        else { std::cerr << "\nte-debug, variable [" << sk2 << "] not found\n"; }
#endif
			}
		    }
		}
	    return oss.str();
	    }

	/**
	 * @brief Fills the content of the block to their destination
	 *
	 * @param crsPage The input template text
	 * @param mData The variables to fill in
	 */
	std::string FillBlocks(std::string const & crsPage, TMapOfBlocks const & mBlocks) const
	    {
	    std::ostringstream oss{};

	    std::regex const re(R"(\{\%\s+block\s+|\{\%\s+endblock\s+\%\})");
	    size_t n{0};
	    for (auto it = std::sregex_token_iterator(crsPage.begin(), crsPage.end(), re, -1); it != std::sregex_token_iterator(); ++it)
		{
		if ( ++n & 1 )
		    {
		    oss << *it;
		    }
		else
		    {
		    std::string const s = *it;
		    size_t      const p = s.find(' ');
		    std::string const k = s.substr(0, p);
		    size_t      const q = s.find('}')+1;
                    if ( mBlocks.find(k) != mBlocks.end() )
                        {
		        oss << mBlocks.at(k);
                        }
                    else
                        {
#ifdef DEBUG
                        std::cerr << "\nte-debug, block [" << k << "] not overwritten, keep: [" << s.substr(q) << "]\n";
#endif
                        oss << s.substr(q);
                        }
		    }
		}
	    return oss.str();
	    }

	/**
	 * @brief Collects the content of all blocks.
	 *
	 * If a template contains block, all other content can be ignores besides the
	 * function 'extends' to find a shell to put the blocks into
	 *
	 * @param page The template
	 */
	TMapOfBlocks GetBlocks(std::string const & crsPage) const
	    {
	    TMapOfBlocks mResult{};

	    std::regex const re(R"(\{\%\s+block\s+|\{\%\s+endblock\s+\%\})");
	    size_t n{0};
	    for (auto it = std::sregex_token_iterator(crsPage.begin(), crsPage.end(), re, -1); it != std::sregex_token_iterator(); ++it)
		{
		if ( n++ & 1 )
		    {
		    std::string const s = *it;
		    size_t      const p = s.find(' ');
		    size_t      const q = s.find('}');
		    mResult[s.substr(0, p)] = s.substr(q+1);
		    }
		}
	    return std::move(mResult);
	    }

	/**
	 * @brief Filters the template by removing sequences surrounded by
	 * 	  surrounded by negative evaluated IF statements like this:
	 *
	 * 	  "{%if data %}{% endif%}" becomes ""
	 *
	 * 	  This enables to remove prefix and suffix of loops
	 */
	std::string GetIfs(std::string const & crsPage, TRenderData const & mData) const
	    {
	    std::ostringstream oss{};

	    std::regex const re(R"(\{\%\s+if\s+|\{\%\s+endif\s+\%\})");
	    size_t n{0};
	    for (auto it = std::sregex_token_iterator(crsPage.begin(), crsPage.end(), re, -1); it != std::sregex_token_iterator(); ++it)
		{
		if ( ++n & 1 )
                    {
                    oss << *it;
                    }
                else
		    {
		    std::string const s = *it;
		    size_t      const p = s.find(' ');
		    size_t      const q = s.find('}');
		    if ( mData.find(s.substr(0, p)) != mData.end() )
                        {
                        oss << s.substr(q+1);
                        }
		    }
		}
	    return std::move(oss.str());
	    }

	/**
	 * @brief Reads the given file into a string
	 *
	 *
	 */
	std::string ReadTemplate(std::string const & crsFilename) const
	    {
	    std::ifstream f(crsFilename);
	    std::string   t((std::istreambuf_iterator<char>(f)),
			     std::istreambuf_iterator<char>());
	    return std::move(t);
	    }

    }; // class CTE

// __TE_H
#endif
