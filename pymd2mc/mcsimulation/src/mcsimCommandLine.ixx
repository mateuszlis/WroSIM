// This code was generated by CLI, a command line interface
// compiler for C++.
//

namespace cli
{
  // unknown_mode
  //
  inline unknown_mode::
  unknown_mode (value v)
  : v_ (v)
  {
  }

  // exception
  //
  inline std::ostream&
  operator<< (std::ostream& os, const exception& e)
  {
    e.print (os);
    return os;
  }

  // unknown_option
  //
  inline unknown_option::
  unknown_option (const std::string& option)
  : option_ (option)
  {
  }

  inline const std::string& unknown_option::
  option () const
  {
    return option_;
  }

  // unknown_argument
  //
  inline unknown_argument::
  unknown_argument (const std::string& argument)
  : argument_ (argument)
  {
  }

  inline const std::string& unknown_argument::
  argument () const
  {
    return argument_;
  }

  // missing_value
  //
  inline missing_value::
  missing_value (const std::string& option)
  : option_ (option)
  {
  }

  inline const std::string& missing_value::
  option () const
  {
    return option_;
  }

  // invalid_value
  //
  inline invalid_value::
  invalid_value (const std::string& option,
                 const std::string& value)
  : option_ (option),  value_ (value)
  {
  }

  inline const std::string& invalid_value::
  option () const
  {
    return option_;
  }

  inline const std::string& invalid_value::
  value () const
  {
    return value_;
  }

  // argv_scanner
  //
  inline argv_scanner::
  argv_scanner (int& argc, char** argv, bool erase)
  : i_ (1), argc_ (argc), argv_ (argv), erase_ (erase)
  {
  }

  inline argv_scanner::
  argv_scanner (int start, int& argc, char** argv, bool erase)
  : i_ (start), argc_ (argc), argv_ (argv), erase_ (erase)
  {
  }

  inline int argv_scanner::
  end () const
  {
    return i_;
  }
}

// options
//

inline const bool& options::
help () const
{
  return this->help_;
}

inline const std::string& options::
sampling () const
{
  return this->sampling_;
}

inline const std::string& options::
f () const
{
  return this->f_;
}

inline const int& options::
first_particles_count () const
{
  return this->first_particles_count_;
}

inline const int& options::
proteins_count () const
{
  return this->proteins_count_;
}

inline const int& options::
latt_row_size () const
{
  return this->latt_row_size_;
}

inline const int& options::
latt_row_count () const
{
  return this->latt_row_count_;
}

inline const std::string& options::
o () const
{
  return this->o_;
}

inline const int& options::
T () const
{
  return this->T_;
}

inline const int& options::
steps () const
{
  return this->steps_;
}

inline const int& options::
output_freq () const
{
  return this->output_freq_;
}

inline const float& options::
omegaAB () const
{
  return this->omegaAB_;
}

inline const float& options::
omegaAC () const
{
  return this->omegaAC_;
}

inline const float& options::
omegaBC () const
{
  return this->omegaBC_;
}

inline const int& options::
eq_steps () const
{
  return this->eq_steps_;
}

inline const bool& options::
no_random_start () const
{
  return this->no_random_start_;
}

inline const bool& options::
enable_calc_msd () const
{
  return this->enable_calc_msd_;
}

