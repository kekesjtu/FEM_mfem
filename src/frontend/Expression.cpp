#include "fem/frontend/Expression.hpp"

#include "muParser.h"

#include <stdexcept>
#include <utility>

namespace fem::frontend
{
Expression::Expression(std::string expression) : expression_(std::move(expression))
{
    try
    {
        size_t processed = 0;
        constant_value_ = std::stod(expression_, &processed);
        is_constant_ = (processed == expression_.size());
    }
    catch (const std::exception &)
    {
        is_constant_ = false;
    }

    if (!is_constant_)
    {
        InitializeParser();
    }
}

Expression::Expression(const Expression &other)
    : expression_(other.expression_),
      is_constant_(other.is_constant_),
      constant_value_(other.constant_value_)
{
    if (!is_constant_)
    {
        InitializeParser();
    }
}

Expression &Expression::operator=(const Expression &other)
{
    if (this == &other)
    {
        return *this;
    }

    expression_ = other.expression_;
    is_constant_ = other.is_constant_;
    constant_value_ = other.constant_value_;
    parser_.reset();
    if (!is_constant_)
    {
        InitializeParser();
    }

    return *this;
}

Expression::Expression(Expression &&other) noexcept
    : expression_(std::move(other.expression_)),
      is_constant_(other.is_constant_),
      constant_value_(other.constant_value_),
      x_(other.x_),
      y_(other.y_),
      z_(other.z_),
      t_(other.t_),
      T_(other.T_),
      parser_(std::move(other.parser_))
{
    if (parser_)
    {
        parser_->DefineVar("x", &x_);
        parser_->DefineVar("y", &y_);
        parser_->DefineVar("z", &z_);
        parser_->DefineVar("t", &t_);
        parser_->DefineVar("T", &T_);
    }
}

Expression &Expression::operator=(Expression &&other) noexcept
{
    if (this == &other)
        return *this;

    expression_ = std::move(other.expression_);
    is_constant_ = other.is_constant_;
    constant_value_ = other.constant_value_;
    x_ = other.x_;
    y_ = other.y_;
    z_ = other.z_;
    t_ = other.t_;
    T_ = other.T_;
    parser_ = std::move(other.parser_);

    // 【修复修复】：同样需要重新绑定地址
    if (parser_)
    {
        parser_->DefineVar("x", &x_);
        parser_->DefineVar("y", &y_);
        parser_->DefineVar("z", &z_);
        parser_->DefineVar("t", &t_);
        parser_->DefineVar("T", &T_);
    }

    return *this;
}

Expression::~Expression() = default;

void Expression::InitializeParser()
{
    parser_ = std::make_unique<mu::Parser>();
    parser_->DefineVar("x", &x_);
    parser_->DefineVar("y", &y_);
    parser_->DefineVar("z", &z_);
    parser_->DefineVar("t", &t_);
    parser_->DefineVar("T", &T_);

    try
    {
        parser_->SetExpr(expression_);
    }
    catch (const mu::Parser::exception_type &e)
    {
        throw std::runtime_error("muparser expression parse failed: " + e.GetMsg());
    }
}

double Expression::Evaluate(const EvalContext &ctx) const
{
    if (is_constant_)
    {
        return constant_value_;
    }

    x_ = ctx.x;
    y_ = ctx.y;
    z_ = ctx.z;
    t_ = ctx.t;
    T_ = ctx.T;

    try
    {
        return parser_->Eval();
    }
    catch (const mu::Parser::exception_type &e)
    {
        throw std::runtime_error("muparser evaluation failed: " + e.GetMsg());
    }
}
}  // namespace fem::frontend
