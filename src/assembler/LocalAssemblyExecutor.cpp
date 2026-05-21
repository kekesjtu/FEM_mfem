#include "assembler/LocalAssemblyExecutor.hpp"

namespace fem::assembler
{
namespace
{
bool MarkerContains(const mfem::Array<int> *marker, int bdr_attr)
{
    if (!marker)
    {
        return true;
    }
    if (bdr_attr < 1 || bdr_attr > marker->Size())
    {
        return false;
    }
    return (*marker)[bdr_attr - 1] != 0;
}

std::vector<std::unique_ptr<mfem::BilinearFormIntegrator>> CreateBilinearIntegrators(
    const std::vector<BilinearIntegratorFactory> &factories)
{
    std::vector<std::unique_ptr<mfem::BilinearFormIntegrator>> integrators;
    integrators.reserve(factories.size());
    for (const auto &factory : factories)
    {
        integrators.push_back(factory());
    }
    return integrators;
}

std::vector<std::unique_ptr<mfem::LinearFormIntegrator>> CreateLinearIntegrators(
    const std::vector<LinearIntegratorFactory> &factories)
{
    std::vector<std::unique_ptr<mfem::LinearFormIntegrator>> integrators;
    integrators.reserve(factories.size());
    for (const auto &factory : factories)
    {
        integrators.push_back(factory());
    }
    return integrators;
}
}  // namespace

std::unique_ptr<mfem::ParBilinearForm> AssembleBilinearFormOpenMP(
    mfem::ParFiniteElementSpace &space,
    const std::vector<BilinearIntegratorFactory> &domain_integrators,
    const std::vector<BoundaryBilinearIntegratorSpec> &boundary_integrators)
{
    auto form = std::make_unique<mfem::ParBilinearForm>(&space);

    const int ne = space.GetNE();
    std::vector<mfem::DenseMatrix> element_matrices(ne);  // 存放ne个单元的单元矩阵

    // C++ 标准对 std::vector<bool> 做了特化.它不是存 bool的数组，而是把每个元素压成一个
    // bit（位压缩存储）。这会带来问题。因此这里使用 std::vector<char> 来存储标志。
    std::vector<char> has_element_matrix(ne, 0);

#pragma omp parallel if (ne > 0)
    {
        auto local_domain_integrators = CreateBilinearIntegrators(domain_integrators);
        mfem::DenseMatrix sum;
        mfem::DenseMatrix tmp;
        mfem::IsoparametricTransformation trans;

#pragma omp for schedule(static)
        for (int i = 0; i < ne; ++i)  // 遍历ne个单元
        {
            if (local_domain_integrators.empty())
            {
                continue;
            }

            space.GetElementTransformation(i, &trans);  // 获取第i个单元的变换规则

            bool has_sum = false;
            sum.SetSize(0);
            const mfem::FiniteElement &fe = *space.GetFE(i);  // 获取第i个单元的有限元信息
            for (auto &integrator : local_domain_integrators)
            {
                // 分别使用local_domain_integrators（积分子的vector）中
                // 各类积分子，组装第i个单元的单元矩阵
                integrator->AssembleElementMatrix(fe, trans, tmp);
                if (tmp.Size() == 0)
                {
                    continue;
                }
                if (!has_sum)
                {
                    sum = tmp;
                    has_sum = true;
                }
                else
                {
                    sum += tmp;
                }
            }

            if (has_sum)
            {
                element_matrices[i] = sum;
                has_element_matrix[i] = 1;
            }
        }
    }

    mfem::Array<int> vdofs;
    mfem::DofTransformation doftrans;
    for (int i = 0; i < ne; ++i)
    {
        if (has_element_matrix[i])
        {
            space.GetElementVDofs(i, vdofs, doftrans);
            doftrans.TransformDual(element_matrices[i]);
            // 将第i个单元的单元矩阵组装到全局矩阵中
            form->AssembleElementMatrix(i, element_matrices[i]);
        }
    }

    if (!boundary_integrators.empty())
    {
        const int nbe = space.GetNBE();
        std::vector<mfem::DenseMatrix> boundary_matrices(nbe);
        std::vector<char> has_boundary_matrix(nbe, 0);

#pragma omp parallel if (nbe > 0)
        {
            std::vector<std::unique_ptr<mfem::BilinearFormIntegrator>> local_boundary_integrators;
            local_boundary_integrators.reserve(boundary_integrators.size());
            for (const auto &spec : boundary_integrators)
            {
                local_boundary_integrators.push_back(spec.factory());
            }

            mfem::DenseMatrix sum;
            mfem::DenseMatrix tmp;
            mfem::IsoparametricTransformation trans;

#pragma omp for schedule(static)
            for (int boundary_element = 0; boundary_element < nbe; ++boundary_element)
            {
                const int bdr_attr = space.GetBdrAttribute(boundary_element);
                space.GetParMesh()->GetBdrElementTransformation(boundary_element, &trans);

                bool has_sum = false;
                sum.SetSize(0);
                const mfem::FiniteElement &be = *space.GetBE(boundary_element);
                for (std::size_t k = 0; k < boundary_integrators.size(); ++k)
                {
                    if (!MarkerContains(boundary_integrators[k].marker, bdr_attr))
                    {
                        continue;
                    }
                    local_boundary_integrators[k]->AssembleElementMatrix(be, trans, tmp);
                    if (tmp.Size() == 0)
                    {
                        continue;
                    }
                    if (!has_sum)
                    {
                        sum = tmp;
                        has_sum = true;
                    }
                    else
                    {
                        sum += tmp;
                    }
                }

                if (has_sum)
                {
                    boundary_matrices[boundary_element] = sum;
                    has_boundary_matrix[boundary_element] = 1;
                }
            }
        }

        for (int boundary_element = 0; boundary_element < nbe; ++boundary_element)
        {
            if (has_boundary_matrix[boundary_element])
            {
                space.GetBdrElementVDofs(boundary_element, vdofs, doftrans);
                doftrans.TransformDual(boundary_matrices[boundary_element]);
                form->AssembleBdrElementMatrix(boundary_element,
                                               boundary_matrices[boundary_element]);
            }
        }
    }

    return form;
}

std::unique_ptr<mfem::ParLinearForm> AssembleLinearFormOpenMP(
    mfem::ParFiniteElementSpace &space,
    const std::vector<LinearIntegratorFactory> &domain_integrators,
    const std::vector<BoundaryLinearIntegratorSpec> &boundary_integrators)
{
    auto form = std::make_unique<mfem::ParLinearForm>(&space);
    *form = 0.0;

    const int ne = space.GetNE();
    std::vector<mfem::Vector> element_vectors(ne);
    std::vector<char> has_element_vector(ne, 0);

#pragma omp parallel if (ne > 0)
    {
        auto local_domain_integrators = CreateLinearIntegrators(domain_integrators);
        mfem::Vector sum;
        mfem::Vector tmp;
        mfem::IsoparametricTransformation trans;

#pragma omp for schedule(static)
        for (int i = 0; i < ne; ++i)
        {
            if (local_domain_integrators.empty())
            {
                continue;
            }

            space.GetElementTransformation(i, &trans);

            bool has_sum = false;
            sum.SetSize(0);
            const mfem::FiniteElement &fe = *space.GetFE(i);
            for (auto &integrator : local_domain_integrators)
            {
                integrator->AssembleRHSElementVect(fe, trans, tmp);
                if (tmp.Size() == 0)
                {
                    continue;
                }
                if (!has_sum)
                {
                    sum = tmp;
                    has_sum = true;
                }
                else
                {
                    sum += tmp;
                }
            }

            if (has_sum)
            {
                element_vectors[i] = sum;
                has_element_vector[i] = 1;
            }
        }
    }

    mfem::Array<int> vdofs;
    mfem::DofTransformation doftrans;
    for (int i = 0; i < ne; ++i)
    {
        if (has_element_vector[i])
        {
            // 这里线性型将单元向量累加到全局向量时，不能像双线性型那样直接直接通过单元索引组装，
            // 而是得先获取第i个单元的vdofs自由度向量才能累加。
            // 纯粹因为mfem::LinearForm的接口设计不对称。
            space.GetElementVDofs(i, vdofs, doftrans);
            doftrans.TransformDual(element_vectors[i]);
            form->AddElementVector(vdofs, element_vectors[i]);
        }
    }

    if (!boundary_integrators.empty())
    {
        const int nbe = space.GetNBE();
        std::vector<mfem::Vector> boundary_vectors(nbe);
        std::vector<char> has_boundary_vector(nbe, 0);

#pragma omp parallel if (nbe > 0)
        {
            std::vector<std::unique_ptr<mfem::LinearFormIntegrator>> local_boundary_integrators;
            local_boundary_integrators.reserve(boundary_integrators.size());
            for (const auto &spec : boundary_integrators)
            {
                local_boundary_integrators.push_back(spec.factory());
            }

            mfem::Vector sum;
            mfem::Vector tmp;
            mfem::IsoparametricTransformation trans;

#pragma omp for schedule(static)
            for (int boundary_element = 0; boundary_element < nbe; ++boundary_element)
            {
                const int bdr_attr = space.GetBdrAttribute(boundary_element);
                space.GetParMesh()->GetBdrElementTransformation(boundary_element, &trans);

                bool has_sum = false;
                sum.SetSize(0);
                const mfem::FiniteElement &be = *space.GetBE(boundary_element);
                for (std::size_t k = 0; k < boundary_integrators.size(); ++k)
                {
                    if (!MarkerContains(boundary_integrators[k].marker, bdr_attr))
                    {
                        continue;
                    }
                    local_boundary_integrators[k]->AssembleRHSElementVect(be, trans, tmp);
                    if (tmp.Size() == 0)
                    {
                        continue;
                    }
                    if (!has_sum)
                    {
                        sum = tmp;
                        has_sum = true;
                    }
                    else
                    {
                        sum += tmp;
                    }
                }

                if (has_sum)
                {
                    boundary_vectors[boundary_element] = sum;
                    has_boundary_vector[boundary_element] = 1;
                }
            }
        }

        mfem::Array<int> boundary_vdofs;
        for (int boundary_element = 0; boundary_element < nbe; ++boundary_element)
        {
            if (has_boundary_vector[boundary_element])
            {
                space.GetBdrElementVDofs(boundary_element, boundary_vdofs, doftrans);
                doftrans.TransformDual(boundary_vectors[boundary_element]);
                form->AddElementVector(boundary_vdofs, boundary_vectors[boundary_element]);
            }
        }
    }

    return form;
}

}  // namespace fem::assembler
