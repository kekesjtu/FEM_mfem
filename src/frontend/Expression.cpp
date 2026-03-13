#include "fem/frontend/Expression.hpp"

#include "muParser.h"

#include <stdexcept>
#include <utility>

namespace fem::frontend
{
ScalarExpression::ScalarExpression(std::string expression,
                                   std::vector<std::string> coupled_variable_names)
    : expression_(std::move(expression)), coupled_variable_names_(std::move(coupled_variable_names))
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

ScalarExpression::ScalarExpression(const ScalarExpression &other)
    : expression_(other.expression_),
      coupled_variable_names_(other.coupled_variable_names_),
      is_constant_(other.is_constant_),
      constant_value_(other.constant_value_)
{
    if (!is_constant_)
    {
        InitializeParser();
    }
}

ScalarExpression &ScalarExpression::operator=(const ScalarExpression &other)
{
    if (this == &other)
    {
        return *this;
    }

    expression_ = other.expression_;
    coupled_variable_names_ = other.coupled_variable_names_;
    is_constant_ = other.is_constant_;
    constant_value_ = other.constant_value_;
    parser_.reset();
    coupled_storage_.clear();

    if (!is_constant_)
    {
        InitializeParser();
    }

    return *this;
}

ScalarExpression::ScalarExpression(ScalarExpression &&other) noexcept
    : expression_(std::move(other.expression_)),
      coupled_variable_names_(std::move(other.coupled_variable_names_)),
      is_constant_(other.is_constant_),
      constant_value_(other.constant_value_),
      x_(other.x_),
      y_(other.y_),
      z_(other.z_),
      t_(other.t_),
      coupled_storage_(std::move(other.coupled_storage_)),
      parser_(std::move(other.parser_))
{
}

ScalarExpression &ScalarExpression::operator=(ScalarExpression &&other) noexcept
{
    if (this == &other)
    {
        return *this;
    }

    expression_ = std::move(other.expression_);
    coupled_variable_names_ = std::move(other.coupled_variable_names_);
    is_constant_ = other.is_constant_;
    constant_value_ = other.constant_value_;
    x_ = other.x_;
    y_ = other.y_;
    z_ = other.z_;
    t_ = other.t_;
    coupled_storage_ = std::move(other.coupled_storage_);
    parser_ = std::move(other.parser_);

    return *this;
}

ScalarExpression::~ScalarExpression() = default;

void ScalarExpression::InitializeParser()
{
    parser_ = std::make_unique<mu::Parser>();
    parser_->DefineVar("x", &x_);
    parser_->DefineVar("y", &y_);
    parser_->DefineVar("z", &z_);
    parser_->DefineVar("t", &t_);

    for (const auto &name : coupled_variable_names_)
    {
        if (name == "x" || name == "y" || name == "z" || name == "t")
        {
            continue;
        }
        auto [it, _] = coupled_storage_.emplace(name, 0.0);
        parser_->DefineVar(name, &it->second);
    }

    try
    {
        parser_->SetExpr(expression_);
    }
    catch (const mu::Parser::exception_type &e)
    {
        throw std::runtime_error("muparser expression setup failed: " + e.GetMsg());
    }
}

double ScalarExpression::Evaluate(const EvalContext &ctx) const
{
    if (is_constant_)
    {
        return constant_value_;
    }

    x_ = ctx.x;
    y_ = ctx.y;
    z_ = ctx.z;
    t_ = ctx.t;

    for (auto &[name, value] : coupled_storage_)
    {
        const auto it = ctx.coupled_fields.find(name);
        value = (it != ctx.coupled_fields.end()) ? it->second : 0.0;
    }

    try
    {
        return parser_->Eval();
    }
    catch (const mu::Parser::exception_type &e)
    {
        throw std::runtime_error("muparser evaluation failed: " + e.GetMsg());
    }
}

const std::string &ScalarExpression::Raw() const
{
    return expression_;
}
}  // namespace fem::frontend
