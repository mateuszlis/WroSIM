// This code was generated by CLI, a command line interface
// compiler for C++.
//

#ifndef ANALYZE_PYRKOVA_COMMAND_LINE_HXX
#define ANALYZE_PYRKOVA_COMMAND_LINE_HXX

#include <iosfwd>
#include <string>
#include <exception>

namespace cli
{
  class unknown_mode
  {
    public:
    enum value
    {
      skip,
      stop,
      fail
    };

    unknown_mode (value v);

    operator value () const 
    {
      return v_;
    }

    private:
    value v_;
  };

  // Exceptions.
  //

  class exception: public std::exception
  {
    public:
    virtual void
    print (std::ostream&) const = 0;
  };

  std::ostream&
  operator<< (std::ostream&, const exception&);

  class unknown_option: public exception
  {
    public:
    virtual
    ~unknown_option () throw ();

    unknown_option (const std::string& option);

    const std::string&
    option () const;

    virtual void
    print (std::ostream&) const;

    virtual const char*
    what () const throw ();

    private:
    std::string option_;
  };

  class unknown_argument: public exception
  {
    public:
    virtual
    ~unknown_argument () throw ();

    unknown_argument (const std::string& argument);

    const std::string&
    argument () const;

    virtual void
    print (std::ostream&) const;

    virtual const char*
    what () const throw ();

    private:
    std::string argument_;
  };

  class missing_value: public exception
  {
    public:
    virtual
    ~missing_value () throw ();

    missing_value (const std::string& option);

    const std::string&
    option () const;

    virtual void
    print (std::ostream&) const;

    virtual const char*
    what () const throw ();

    private:
    std::string option_;
  };

  class invalid_value: public exception
  {
    public:
    virtual
    ~invalid_value () throw ();

    invalid_value (const std::string& option,
                   const std::string& value);

    const std::string&
    option () const;

    const std::string&
    value () const;

    virtual void
    print (std::ostream&) const;

    virtual const char*
    what () const throw ();

    private:
    std::string option_;
    std::string value_;
  };

  class eos_reached: public exception
  {
    public:
    virtual void
    print (std::ostream&) const;

    virtual const char*
    what () const throw ();
  };

  class scanner
  {
    public:
    virtual
    ~scanner ();

    virtual bool
    more () = 0;

    virtual const char*
    peek () = 0;

    virtual const char*
    next () = 0;

    virtual void
    skip () = 0;
  };

  class argv_scanner: public scanner
  {
    public:
    argv_scanner (int& argc, char** argv, bool erase = false);
    argv_scanner (int start, int& argc, char** argv, bool erase = false);

    int
    end () const;

    virtual bool
    more ();

    virtual const char*
    peek ();

    virtual const char*
    next ();

    virtual void
    skip ();

    private:
    int i_;
    int& argc_;
    char** argv_;
    bool erase_;
  };
}

#include <string>

class options
{
  public:

  options (int& argc,
           char** argv,
           bool erase = false,
           ::cli::unknown_mode option = ::cli::unknown_mode::fail,
           ::cli::unknown_mode argument = ::cli::unknown_mode::stop);

  options (int start,
           int& argc,
           char** argv,
           bool erase = false,
           ::cli::unknown_mode option = ::cli::unknown_mode::fail,
           ::cli::unknown_mode argument = ::cli::unknown_mode::stop);

  options (int& argc,
           char** argv,
           int& end,
           bool erase = false,
           ::cli::unknown_mode option = ::cli::unknown_mode::fail,
           ::cli::unknown_mode argument = ::cli::unknown_mode::stop);

  options (int start,
           int& argc,
           char** argv,
           int& end,
           bool erase = false,
           ::cli::unknown_mode option = ::cli::unknown_mode::fail,
           ::cli::unknown_mode argument = ::cli::unknown_mode::stop);

  options (::cli::scanner&,
           ::cli::unknown_mode option = ::cli::unknown_mode::fail,
           ::cli::unknown_mode argument = ::cli::unknown_mode::stop);

  // Option accessors.
  //

  const bool&
  help () const;

  const std::string&
  f () const;

  const std::string&
  o () const;

  const int&
  t () const;

  const float&
  d () const;

  // Print usage information.
  //
  static void
  print_usage (::std::ostream&);

  private:
  void
  _parse (::cli::scanner&,
          ::cli::unknown_mode option,
          ::cli::unknown_mode argument);

  public:
  bool help_;
  std::string f_;
  std::string o_;
  int t_;
  float d_;
};

#include "analyzePyrkovaCommandLine.ixx"

#endif // ANALYZE_PYRKOVA_COMMAND_LINE_HXX
