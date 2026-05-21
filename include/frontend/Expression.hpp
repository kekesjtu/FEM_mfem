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

    Expression(const Expression &other) = default;
    Expression &operator=(const Expression &other) = default;
    Expression(Expression &&other) noexcept = default;
    Expression &operator=(Expression &&other) noexcept = default;
    ~Expression() = default;

    double Evaluate(const EvalContext &ctx) const;

    /// Check if the expression depends on a given variable.
    /// For constant expressions, always returns false.
    /// For non-constant expressions, evaluates at two different values and checks if result
    /// changes.
    bool UsesVariable(const std::string &var_name) const;

    const std::string &GetExpressionString() const
    {
        return expression_;
    }

    bool IsConstant() const
    {
        return is_constant_;
    }

    double ConstantValue() const
    {
        return constant_value_;
    }

  private:
    std::string expression_;
    bool is_constant_ = false;
    double constant_value_ = 0.0;
};

class ExpressionRuntime
{
  public:
    explicit ExpressionRuntime(const Expression &expression);
    ExpressionRuntime(const ExpressionRuntime &other) = delete;
    ExpressionRuntime &operator=(const ExpressionRuntime &other) = delete;
    ExpressionRuntime(ExpressionRuntime &&other) noexcept;
    ExpressionRuntime &operator=(ExpressionRuntime &&other) noexcept;
    ~ExpressionRuntime();

    double Evaluate(const EvalContext &ctx);

  private:
    void InitializeParser();
    void RebindVariables();

    std::string expression_;
    bool is_constant_ = false;
    double constant_value_ = 0.0;

    EvalContext eval_context_;
    std::unique_ptr<mu::Parser> parser_;
};
}  // namespace fem::frontend
