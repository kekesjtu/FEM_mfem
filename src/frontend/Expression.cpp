#include "frontend/Expression.hpp"

#include "muParser.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fem::frontend
{
namespace
{
bool TryParseConstant(const std::string &expression, double &value)
{
    try
    {
        size_t processed = 0;
        value = std::stod(expression, &processed);
        return processed == expression.size();
    }
    catch (const std::exception &)
    {
        return false;
    }
}
}  // namespace

Expression::Expression(std::string expression) : expression_(std::move(expression))
{
    is_constant_ = TryParseConstant(expression_, constant_value_);
}

double Expression::Evaluate(const EvalContext &ctx) const
{
    ExpressionRuntime runtime(*this);
    return runtime.Evaluate(ctx);
}

bool Expression::UsesVariable(const std::string &var_name) const
{
    if (is_constant_)
    {
        return false;
    }

    auto make_ctx = [&](double val) -> EvalContext
    {
        EvalContext ctx{1.0, 2.0, 3.0, 0.0, 293.15};
        if (var_name == "x")
            ctx.x = val;
        else if (var_name == "y")
            ctx.y = val;
        else if (var_name == "z")
            ctx.z = val;
        else if (var_name == "t")
            ctx.t = val;
        else if (var_name == "T")
            ctx.T = val;
        return ctx;
    };

    ExpressionRuntime runtime(*this);
    const double v1 = runtime.Evaluate(make_ctx(100.0));
    const double v2 = runtime.Evaluate(make_ctx(200.0));
    return v1 != v2;
}

ExpressionRuntime::ExpressionRuntime(const Expression &expression)
    : expression_(expression.GetExpressionString()),
      is_constant_(expression.IsConstant()),
      constant_value_(expression.ConstantValue())
{
    if (!is_constant_)
    {
        InitializeParser();
    }
}

ExpressionRuntime::ExpressionRuntime(ExpressionRuntime &&other) noexcept
    : expression_(std::move(other.expression_)),
      is_constant_(other.is_constant_),
      constant_value_(other.constant_value_),
      eval_context_(other.eval_context_),
      parser_(std::move(other.parser_))
{
    RebindVariables();
}

ExpressionRuntime &ExpressionRuntime::operator=(ExpressionRuntime &&other) noexcept
{
    if (this == &other)
    {
        return *this;
    }

    expression_ = std::move(other.expression_);
    is_constant_ = other.is_constant_;
    constant_value_ = other.constant_value_;
    eval_context_ = other.eval_context_;
    parser_ = std::move(other.parser_);
    RebindVariables();
    return *this;
}

ExpressionRuntime::~ExpressionRuntime() = default;

void ExpressionRuntime::InitializeParser()
{
    parser_ = std::make_unique<mu::Parser>();
    RebindVariables();
    parser_->DefineConst("pi", M_PI);

    try
    {
        parser_->SetExpr(expression_);
    }
    catch (const mu::Parser::exception_type &e)
    {
        throw std::runtime_error("muparser expression parse failed: " + e.GetMsg());
    }
}

void ExpressionRuntime::RebindVariables()
{
    if (!parser_)
    {
        return;
    }

    parser_->DefineVar("x", &eval_context_.x);
    parser_->DefineVar("y", &eval_context_.y);
    parser_->DefineVar("z", &eval_context_.z);
    parser_->DefineVar("t", &eval_context_.t);
    parser_->DefineVar("T", &eval_context_.T);
}

double ExpressionRuntime::Evaluate(const EvalContext &ctx)
{
    if (is_constant_)
    {
        return constant_value_;
    }

    eval_context_ = ctx;

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
