#pragma once

#include <memory>
#include <string>

namespace mu
{
class Parser;
}

namespace fem::frontend
{
struct EvalContext
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double t = 0.0;
    double T = 0.0;
};

class Expression
{
  public:
    explicit Expression(std::string expression);

    Expression(const Expression &other);
    Expression &operator=(const Expression &other);
    Expression(Expression &&other) noexcept;
    Expression &operator=(Expression &&other) noexcept;
    ~Expression();

    double Evaluate(const EvalContext &ctx) const;

  private:
    void InitializeParser();

    std::string expression_;
    bool is_constant_ = false;
    double constant_value_ = 0.0;

    mutable double x_ = 0.0;
    mutable double y_ = 0.0;
    mutable double z_ = 0.0;
    mutable double t_ = 0.0;
    mutable double T_ = 0.0;
    std::unique_ptr<mu::Parser> parser_;
};
}  // namespace fem::frontend
