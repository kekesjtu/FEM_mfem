#include "fem/mapping/CoefficientFactory.hpp"

namespace fem::mapping
{
std::unique_ptr<mfem::Coefficient> CoefficientFactory::CreateScalar(
    const fem::frontend::ScalarExpression &expr)
{
    return std::make_unique<mfem::FunctionCoefficient>(
        [expr](const mfem::Vector &x) -> double
        {
            fem::frontend::EvalContext ctx;
            ctx.x = x.Size() > 0 ? x(0) : 0.0;
            ctx.y = x.Size() > 1 ? x(1) : 0.0;
            ctx.z = x.Size() > 2 ? x(2) : 0.0;
            return expr.Evaluate(ctx);
        });
}
}  // namespace fem::mapping
