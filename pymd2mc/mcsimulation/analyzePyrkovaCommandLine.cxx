// This code was generated by CLI, a command line interface
// compiler for C++.
//

#include "analyzePyrkovaCommandLine.hxx"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <ostream>
#include <sstream>

namespace cli
{
  // unknown_option
  //
  unknown_option::
  ~unknown_option () throw ()
  {
  }

  void unknown_option::
  print (std::ostream& os) const
  {
    os << "unknown option '" << option () << "'";
  }

  const char* unknown_option::
  what () const throw ()
  {
    return "unknown option";
  }

  // unknown_argument
  //
  unknown_argument::
  ~unknown_argument () throw ()
  {
  }

  void unknown_argument::
  print (std::ostream& os) const
  {
    os << "unknown argument '" << argument () << "'";
  }

  const char* unknown_argument::
  what () const throw ()
  {
    return "unknown argument";
  }

  // missing_value
  //
  missing_value::
  ~missing_value () throw ()
  {
  }

  void missing_value::
  print (std::ostream& os) const
  {
    os << "missing value for option '" << option () << "'";
  }

  const char* missing_value::
  what () const throw ()
  {
    return "missing option value";
  }

  // invalid_value
  //
  invalid_value::
  ~invalid_value () throw ()
  {
  }

  void invalid_value::
  print (std::ostream& os) const
  {
    os << "invalid value '" << value () << "' for option '"
       << option () << "'";
  }

  const char* invalid_value::
  what () const throw ()
  {
    return "invalid option value";
  }

  // eos_reached
  //
  void eos_reached::
  print (std::ostream& os) const
  {
    os << what ();
  }

  const char* eos_reached::
  what () const throw ()
  {
    return "end of argument stream reached";
  }

  // scanner
  //
  scanner::
  ~scanner ()
  {
  }

  // argv_scanner
  //
  bool argv_scanner::
  more ()
  {
    return i_ < argc_;
  }

  const char* argv_scanner::
  peek ()
  {
    if (i_ < argc_)
      return argv_[i_];
    else
      throw eos_reached ();
  }

  const char* argv_scanner::
  next ()
  {
    if (i_ < argc_)
    {
      const char* r (argv_[i_]);

      if (erase_)
      {
        for (int i (i_ + 1); i < argc_; ++i)
          argv_[i - 1] = argv_[i];

        --argc_;
        argv_[argc_] = 0;
      }
      else
        ++i_;

      return r;
    }
    else
      throw eos_reached ();
  }

  void argv_scanner::
  skip ()
  {
    if (i_ < argc_)
      ++i_;
    else
      throw eos_reached ();
  }

  template <typename X>
  struct parser
  {
    static void
    parse (X& x, scanner& s)
    {
      const char* o (s.next ());

      if (s.more ())
      {
        const char* v (s.next ());
        std::istringstream is (v);
        if (!(is >> x && is.eof ()))
          throw invalid_value (o, v);
      }
      else
        throw missing_value (o);
    }
  };

  template <>
  struct parser<bool>
  {
    static void
    parse (bool& x, scanner& s)
    {
      s.next ();
      x = true;
    }
  };

  template <>
  struct parser<std::string>
  {
    static void
    parse (std::string& x, scanner& s)
    {
      const char* o (s.next ());

      if (s.more ())
        x = s.next ();
      else
        throw missing_value (o);
    }
  };

  template <typename X>
  struct parser<std::vector<X> >
  {
    static void
    parse (std::vector<X>& c, scanner& s)
    {
      X x;
      parser<X>::parse (x, s);
      c.push_back (x);
    }
  };

  template <typename X>
  struct parser<std::set<X> >
  {
    static void
    parse (std::set<X>& c, scanner& s)
    {
      X x;
      parser<X>::parse (x, s);
      c.insert (x);
    }
  };

  template <typename K, typename V>
  struct parser<std::map<K, V> >
  {
    static void
    parse (std::map<K, V>& m, scanner& s)
    {
      const char* o (s.next ());

      if (s.more ())
      {
        std::string ov (s.next ());
        std::string::size_type p = ov.find ('=');

        if (p == std::string::npos)
        {
          K k = K ();

          if (!ov.empty ())
          {
            std::istringstream ks (ov);

            if (!(ks >> k && ks.eof ()))
              throw invalid_value (o, ov);
          }

          m[k] = V ();
        }
        else
        {
          K k = K ();
          V v = V ();
          std::string kstr (ov, 0, p);
          std::string vstr (ov, p + 1);

          if (!kstr.empty ())
          {
            std::istringstream ks (kstr);

            if (!(ks >> k && ks.eof ()))
              throw invalid_value (o, ov);
          }

          if (!vstr.empty ())
          {
            std::istringstream vs (vstr);

            if (!(vs >> v && vs.eof ()))
              throw invalid_value (o, ov);
          }

          m[k] = v;
        }
      }
      else
        throw missing_value (o);
    }
  };

  template <typename X, typename T, T X::*P>
  void
  thunk (X& x, scanner& s)
  {
    parser<T>::parse (x.*P, s);
  }
}

#include <map>
#include <cstring>

// options
//

options::
options (int& argc,
         char** argv,
         bool erase,
         ::cli::unknown_mode opt,
         ::cli::unknown_mode arg)
: help_ (),
  f_ (),
  t_ (1),
  d_ (1)
{
  ::cli::argv_scanner s (argc, argv, erase);
  _parse (s, opt, arg);
}

options::
options (int start,
         int& argc,
         char** argv,
         bool erase,
         ::cli::unknown_mode opt,
         ::cli::unknown_mode arg)
: help_ (),
  f_ (),
  t_ (1),
  d_ (1)
{
  ::cli::argv_scanner s (start, argc, argv, erase);
  _parse (s, opt, arg);
}

options::
options (int& argc,
         char** argv,
         int& end,
         bool erase,
         ::cli::unknown_mode opt,
         ::cli::unknown_mode arg)
: help_ (),
  f_ (),
  t_ (1),
  d_ (1)
{
  ::cli::argv_scanner s (argc, argv, erase);
  _parse (s, opt, arg);
  end = s.end ();
}

options::
options (int start,
         int& argc,
         char** argv,
         int& end,
         bool erase,
         ::cli::unknown_mode opt,
         ::cli::unknown_mode arg)
: help_ (),
  f_ (),
  t_ (1),
  d_ (1)
{
  ::cli::argv_scanner s (start, argc, argv, erase);
  _parse (s, opt, arg);
  end = s.end ();
}

options::
options (::cli::scanner& s,
         ::cli::unknown_mode opt,
         ::cli::unknown_mode arg)
: help_ (),
  f_ (),
  t_ (1),
  d_ (1)
{
  _parse (s, opt, arg);
}

void options::
print_usage (::std::ostream& os)
{
  os << "--help        Display this message" << ::std::endl;

  os << "-f <filename> Input file [gro]" << ::std::endl;

  os << "-t <num>      Value describes number of frames that are the threshold for" << ::std::endl
     << ::std::endl
     << "              clusters identification" << ::std::endl;

  os << "-d <float>    Distance to identify a clusterized atoms" << ::std::endl;
}

typedef
std::map<std::string, void (*) (options&, ::cli::scanner&)>
_cli_options_map;

static _cli_options_map _cli_options_map_;

struct _cli_options_map_init
{
  _cli_options_map_init ()
  {
    _cli_options_map_["--help"] = 
    &::cli::thunk< options, bool, &options::help_ >;
    _cli_options_map_["-f"] = 
    &::cli::thunk< options, std::string, &options::f_ >;
    _cli_options_map_["-t"] = 
    &::cli::thunk< options, int, &options::t_ >;
    _cli_options_map_["-d"] = 
    &::cli::thunk< options, float, &options::d_ >;
  }
} _cli_options_map_init_;

void options::
_parse (::cli::scanner& s,
        ::cli::unknown_mode opt_mode,
        ::cli::unknown_mode arg_mode)
{
  bool opt = true;

  while (s.more ())
  {
    const char* o = s.peek ();

    if (std::strcmp (o, "--") == 0)
    {
      s.skip ();
      opt = false;
      continue;
    }

    _cli_options_map::const_iterator i (
      opt ? _cli_options_map_.find (o) : _cli_options_map_.end ());

    if (i != _cli_options_map_.end ())
    {
      (*(i->second)) (*this, s);
    }
    else if (opt && std::strncmp (o, "-", 1) == 0 && o[1] != '\0')
    {
      switch (opt_mode)
      {
        case ::cli::unknown_mode::skip:
        {
          s.skip ();
          continue;
        }
        case ::cli::unknown_mode::stop:
        {
          break;
        }
        case ::cli::unknown_mode::fail:
        {
          throw ::cli::unknown_option (o);
        }
      }

      break;
    }
    else
    {
      switch (arg_mode)
      {
        case ::cli::unknown_mode::skip:
        {
          s.skip ();
          continue;
        }
        case ::cli::unknown_mode::stop:
        {
          break;
        }
        case ::cli::unknown_mode::fail:
        {
          throw ::cli::unknown_argument (o);
        }
      }

      break;
    }
  }
}

