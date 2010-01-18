class header {
public:
// Docstrings might be of use...
#if 0
  /// Construct a header with an empty type code and no fields
  header() { }

  /// Construct a header by splitting up a tab-separated text string
  explicit header(const std::string& line) { assign(line); }

  /// @brief Allocate a new particular kind of header by splitting up
  /// a tab-separated text string
  static header* new_header(const std::string& line);

  /// Destroy this header object
  virtual ~header() { }

  /// Returns the field with the specified @a tag, or @c end() if not found.
  iterator find(const char* tag);
#endif

  /// Assign to a header by splitting up a tab-separated text string
  header& assign(const std::string& line);

  /// @brief Text of the field with the given @a tag, or @a default_value
  /// if there is no such field
  std::string field(const char* tag,
		    const std::string& default_value = std::string()) const;

  /// @brief The field with the given @a tag as an integer, or @a default_value
  /// if there is no such field
  int field_int(const char* tag, int default_value = 0) const;


  void set_type(const std::string& type) { type_ = type; sync(); }

  void set_field(const char* tag, const std::string& value);
  void set_field(const char* tag, int value);

  void set_field(iterator pos, const std::string& value);
  void set_field(iterator pos, int value);
};
