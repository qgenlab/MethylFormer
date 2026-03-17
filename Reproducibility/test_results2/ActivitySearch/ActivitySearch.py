import sqlite3
import pandas as pd

class ActivitySearch:
    def __init__(self, db_path = "methylation_benchmark.db"):
        self._conn = sqlite3.connect(db_path)
        self._QUERY_TISSUES_LIST = "SELECT DISTINCT tissue FROM samples ORDER BY tissue"
        self._QUERY_TISSUES_LIST_INFO = "SELECT * FROM samples ORDER BY tissue, sex, age"
        self._QUERY_GENES_LIST = "SELECT * FROM genes ORDER BY gene_name"
        self._QUERY_ACTIVE_FEATURES = """SELECT DISTINCT f.chrom, f.start, f.end, f.feature_type, f.feature_id, g.gene_name, g.gene_id 
        FROM activity a 
        JOIN samples s ON a.sample_id = s.sample_id 
        JOIN features f ON a.feature_id = f.feature_id 
        LEFT JOIN genes g ON f.gene_id = g.gene_id WHERE s.tissue = ?  AND a.activity_state = 'ACTIVE'
        """
        self._QUERY_DIFFERENTIAL_FEATURES = """ SELECT f.feature_id, f.feature_type, f.chrom, f.start, f.end, g.gene_name,
        s1.tissue as tissue_A, a1.activity_state as state_A,
        s2.tissue as tissue_B, a2.activity_state as state_B
        FROM features f
        JOIN activity a1 ON f.feature_id = a1.feature_id
        JOIN samples s1 ON a1.sample_id = s1.sample_id
        JOIN activity a2 ON f.feature_id = a2.feature_id
        JOIN samples s2 ON a2.sample_id = s2.sample_id
        LEFT JOIN genes g ON f.gene_id = g.gene_id
        WHERE s1.tissue = ? AND s2.tissue = ?
        AND a1.activity_state != a2.activity_state
        """
        self._QUERY_INACTIVE_FEATURES = """
        SELECT DISTINCT f.chrom, f.start, f.end, f.feature_type, f.feature_id, g.gene_name, g.gene_id
        FROM activity a
        JOIN samples s ON a.sample_id = s.sample_id
        JOIN features f ON a.feature_id = f.feature_id
        LEFT JOIN genes g ON f.gene_id = g.gene_id
        WHERE s.tissue = ? AND a.activity_state = 'INACTIVE'
        """
        print(f"Successfully connected to: {db_path}")
        print("Database Schema:")
        tables = pd.read_sql("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name", self._conn)
        for table_name in tables['name']:
            print(f"\n[Table: {table_name}]")
            cols = pd.read_sql(f"PRAGMA table_info({table_name})", self._conn)
            for _, row in cols.iterrows():
                pk_marker = " [PK]" if row['pk'] > 0 else ""
                print(f"  - {row['name']} ({row['type']}){pk_marker}")
        print("\n" + "-"*30 + "\n")
    def list_tissues(self):
        """Return unique list of available tissues."""
        current_query = self._QUERY_TISSUES_LIST
        df = pd.read_sql(current_query, self._conn)
        return df['tissue'].tolist()
    def list_tissues_info(self):
        """Return detailed info (tissue, sex, age) for all samples."""
        current_query = self._QUERY_TISSUES_LIST_INFO
        df = pd.read_sql(current_query, self._conn)
        return df
    def list_genes(self):
        """List available genes (id and name)."""
        current_query = self._QUERY_GENES_LIST
        df = pd.read_sql(current_query, self._conn)
        return df
    def get_active_features(self, tissue, sex=None, age=None):
        """
        Get active regions for a specific tissue.
        - If sex/age are None: Ignores them (returns active regions from ANY donor of that tissue).
        - If sex/age provided: Filters for that specific donor.
        """
        current_query = self._QUERY_ACTIVE_FEATURES
        params = [tissue]
        if sex:
            current_query += " AND s.sex = ?"
            params.append(sex)
        if age:
            current_query += " AND s.age = ?"
            params.append(age)
        df = pd.read_sql(current_query, self._conn, params=params)
        return df
    def get_differential_features(self, tissue_A, tissue_B, 
                                  sex_A=None, age_A=None, 
                                  sex_B=None, age_B=None):
        """
        Find features where activity differs between Tissue A and Tissue B.
        - Returns rows where State(A) != State(B).
        - Can filter specifically by sex/age for either A or B.
        """
        current_query = self._QUERY_DIFFERENTIAL_FEATURES
        params = [tissue_A, tissue_B]
        if sex_A:
            current_query += " AND s1.sex = ?"
            params.append(sex_A)
        if age_A:
            current_query += " AND s1.age = ?"
            params.append(age_A)
        if sex_B:
            current_query += " AND s2.sex = ?"
            params.append(sex_B)
        if age_B:
            current_query += " AND s2.age = ?"
            params.append(age_B)
        df = pd.read_sql(current_query, self._conn, params=params)
        return df
    def get_inactive_features(self, tissue, sex=None, age=None):
        """
        Get INACTIVE regions for a specific tissue.
        """
        current_query = self._QUERY_INACTIVE_FEATURES
        params = [tissue]
        if sex:
            current_query += " AND s.sex = ?"
            params.append(sex)
        if age:
            current_query += " AND s.age = ?"
            params.append(age)
        df = pd.read_sql(current_query, self._conn, params=params)
        return df
    def run_custom_query(self, query, params=None):
        """
        Execute a raw SQL query safely.
        Args:
            query (str): The SQL query string. Use '?' for placeholders.
            params (list/tuple): The values for the placeholders (optional).
        Returns:
            pd.DataFrame: Results of the query.
        """
        if params is None:
            params = []
        try:
            df = pd.read_sql(query, self._conn, params=params)
            return df
        except Exception as e:
            print(f"Query failed: {e}")
            return None
    def query_logic(self, formula):
        """
        Parses a logical string to find features.
        Supports:
          - '|' for OR (Union)
          - '&' for AND (Intersection)
          - '-' for NOT (Difference)
          - Parentheses () for grouping

        Example: "liver & lung - kidney"
        (Active in Liver AND Lung, but NOT Kidney)
        """
        import re
        available_tissues = self.list_tissues()
        available_tissues.sort(key=len, reverse=True)
        token_map = {}
        safe_formula = formula
        for i, tissue in enumerate(available_tissues):
            pattern = re.escape(tissue)
            if re.search(pattern, safe_formula, re.IGNORECASE):
                token = f"##TISSUE_{i}##"
                token_map[token] = tissue
                safe_formula = re.sub(pattern, token, safe_formula, flags=re.IGNORECASE)
        safe_formula = safe_formula.replace("|", " UNION ")
        safe_formula = safe_formula.replace("&", " INTERSECT ")
        safe_formula = safe_formula.replace("-", " EXCEPT ")
        base_query = """
        SELECT feature_id FROM activity a 
        JOIN samples s ON a.sample_id = s.sample_id 
        WHERE s.tissue = '{}' AND a.activity_state = 'ACTIVE'
        """
        for token, tissue_name in token_map.items():
            sql_subquery = base_query.format(tissue_name)
            safe_formula = safe_formula.replace(token, sql_subquery)
        final_query = f"""
        WITH logic_result AS (
            {safe_formula}
        )
        SELECT f.chrom, f.start, f.end, f.feature_type, f.feature_id, g.gene_name
        FROM logic_result lr
        JOIN features f ON lr.feature_id = f.feature_id
        LEFT JOIN genes g ON f.gene_id = g.gene_id
        """
        print(f"Executing Logic: {formula}")
        try:
            return pd.read_sql(final_query, self._conn)
        except Exception as e:
            print(f"Logic Error: {e}")
            return None
    def close(self):
        """Helper to close connection when done"""
        self._conn.close()

