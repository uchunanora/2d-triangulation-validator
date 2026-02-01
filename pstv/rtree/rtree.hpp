#pragma once

#include "rtree_fwd.hpp"
#include "rtree_types.hpp"

#include <vector>
#include <algorithm>
#include <iterator>

namespace lib_rtree
{
    /**
     * R-Treeラッパークラス
     * Boost.GeometryのR-Treeをラップして、線分の管理を簡単にする
     */
    template <typename PointType = point_t, typename BoxType = box_t, typename ValueType = rtree_value_t>
    class rtree
    {
    public:
        using point_type = PointType;
        using box_type = BoxType;
        using value_type = ValueType;
        using segment_type = std::pair<std::vector<double>, std::vector<double>>;
        using size_type = std::size_t;

    private:
        // Boost.GeometryのR-Tree（quadratic分割、最大16要素）
        bgi::rtree<value_type, bgi::quadratic<16>> rtree_impl_;

    public:
        /**
         * デフォルトコンストラクタ
         */
        rtree() = default;

        /**
         * デストラクタ
         */
        ~rtree() = default;

        /**
         * コピーコンストラクタ
         */
        rtree(const rtree& other) = default;

        /**
         * ムーブコンストラクタ
         */
        rtree(rtree&& other) noexcept = default;

        /**
         * コピー代入演算子
         */
        rtree& operator=(const rtree& other) = default;

        /**
         * ムーブ代入演算子
         */
        rtree& operator=(rtree&& other) noexcept = default;

        /**
         * 線分を挿入
         * @param p1 線分の端点1
         * @param p2 線分の端点2
         */
        void insert(const std::vector<double>& p1, const std::vector<double>& p2)
        {
            // bounding boxを作成
            box_type box(
                point_type(std::min(p1[0], p2[0]), std::min(p1[1], p2[1])),
                point_type(std::max(p1[0], p2[0]), std::max(p1[1], p2[1]))
            );
            rtree_impl_.insert(std::make_pair(box, std::make_pair(p1, p2)));
        }

        /**
         * 線分を削除
         * @param p1 線分の端点1
         * @param p2 線分の端点2
         */
        void remove(const std::vector<double>& p1, const std::vector<double>& p2)
        {
            // クエリ用のbounding boxを作成
            box_type query_box(
                point_type(std::min(p1[0], p2[0]), std::min(p1[1], p2[1])),
                point_type(std::max(p1[0], p2[0]), std::max(p1[1], p2[1]))
            );
            
            // 交差する可能性のある線分を検索
            std::vector<value_type> results;
            rtree_impl_.query(bgi::intersects(query_box), std::back_inserter(results));
            
            // 一致する線分を削除
            for (const auto& result : results)
            {
                const auto& seg = result.second;
                if (segments_equal(p1, p2, seg.first, seg.second))
                {
                    rtree_impl_.remove(result);
                    break; // 1つの線分だけ削除
                }
            }
        }

        /**
         * 指定されたbounding boxと交差する線分を検索
         * @param p1 クエリ領域の端点1
         * @param p2 クエリ領域の端点2
         * @return 交差する可能性のある線分のリスト
         */
        std::vector<segment_type> query(const std::vector<double>& p1, const std::vector<double>& p2) const
        {
            std::vector<segment_type> candidates;
            
            // クエリ用のbounding boxを作成
            box_type query_box(
                point_type(std::min(p1[0], p2[0]), std::min(p1[1], p2[1])),
                point_type(std::max(p1[0], p2[0]), std::max(p1[1], p2[1]))
            );
            
            // R-treeから交差する可能性のある線分を検索
            std::vector<value_type> results;
            rtree_impl_.query(bgi::intersects(query_box), std::back_inserter(results));
            
            // 結果を変換
            for (const auto& result : results)
            {
                candidates.push_back(result.second);
            }
            
            return candidates;
        }

        /**
         * すべての要素をクリア
         */
        void clear()
        {
            rtree_impl_.clear();
        }

        /**
         * R-treeのサイズを取得
         * @return 格納されている線分の数
         */
        size_type size() const
        {
            return rtree_impl_.size();
        }

        /**
         * R-treeが空かどうかを確認
         * @return 空の場合はtrue
         */
        bool empty() const
        {
            return rtree_impl_.empty();
        }

    private:
        /**
         * 2つの線分が等しいかどうかを判定
         * @param p1 線分1の端点1
         * @param p2 線分1の端点2
         * @param q1 線分2の端点1
         * @param q2 線分2の端点2
         * @return 等しい場合はtrue
         */
        static bool segments_equal(
            const std::vector<double>& p1,
            const std::vector<double>& p2,
            const std::vector<double>& q1,
            const std::vector<double>& q2)
        {
            if (p1 == q1 && p2 == q2)
            {
                return true;
            }
            if (p1 == q2 && p2 == q1)
            {
                return true;
            }
            return false;
        }
    };

    /**
     * デフォルトのR-Tree型
     */
    using rtree_t = rtree<point_t, box_t, rtree_value_t>;
}

