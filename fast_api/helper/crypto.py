from passlib.context import CryptContext

pwd_context = CryptContext(
    schemes=["sha256_crypt", "ldap_salted_md5"],
    sha256_crypt__default_rounds=91234,
    ldap_salted_md5__salt_size=16,
)


def verify_password(plain_password, hashed_password):
    return pwd_context.verify(plain_password, hashed_password)


def get_password_hash(password):
    return pwd_context.hash(password)


if __name__ == "__main__":
    import sys

    password = sys.argv[1]

    hash = get_password_hash(password)

    print(hash)
