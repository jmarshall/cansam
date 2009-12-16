class outer {
public:
  class inner {
  public:
    void do_foo(outer& o);
  };

private:
  void foo() { }
};

void outer::inner::do_foo(outer& o) {
  o.foo();
}
