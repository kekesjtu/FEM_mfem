#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

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
    std::unordered_map<std::string, double> coupled_fields;
};

class ScalarExpression
{
  public:
    explicit ScalarExpression(std::string expression,
                              std::vector<std::string> coupled_variable_names = {});

    ScalarExpression(const ScalarExpression &other);
    ScalarExpression &operator=(const ScalarExpression &other);
    ScalarExpression(ScalarExpression &&other) noexcept;
    ScalarExpression &operator=(ScalarExpression &&other) noexcept;
    ~ScalarExpression();

    double Evaluate(const EvalContext &ctx) const;
    const std::string &Raw() const;

  private:
    void InitializeParser();

    std::string expression_;
    std::vector<std::string> coupled_variable_names_;
    bool is_constant_ = false;
    double constant_value_ = 0.0;

    mutable double x_ = 0.0;
    mutable double y_ = 0.0;
    mutable double z_ = 0.0;
    mutable double t_ = 0.0;
    mutable std::unordered_map<std::string, double> coupled_storage_;
    std::unique_ptr<mu::Parser> parser_;
};
}  // namespace fem::frontend
